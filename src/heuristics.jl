include("utils.jl")

mutable struct HeurVertex
    is_s::Bool
    is_t::Bool
    p::Int64
    ph::Int64
    weight::Float64
end

mutable struct HeurEdge
    d::Float64
    big_d::Float64
    weight::Float64
    src_node::Ref{HeurVertex}
    node_multiplier::Ref{Float64}
end

HeurParams = @NamedTuple begin
    m_start::Float64
    m_smoother::Function
    v_smoother::Function
    e_smoother::Function
end

function baseSmoother(n::Real)::Real
    return n / (n + 1)
end

default_heur_params = (
    m_start=1.0,
    m_smoother=baseSmoother,
    v_smoother=n -> 9 + n / 10 |> baseSmoother,
    e_smoother=n -> 9 + n / 10 |> baseSmoother,
)

function heurWeight(e::HeurEdge)::Float64
    return e.weight + e.node_multiplier[] * e.src_node[].weight
end

# /!\ Prends les s, t du format Graphs.jl, renvoie un path format Graphs.jl (et pas MetaGraphsNext)
function mkPath(state::Graphs.DijkstraState{Float64,Int64}; s::Int64, t::Int64)::Vector{Int64}
    res::Vector = [t]
    while t != s
        t = state.parents[t]
        pushfirst!(res, t)
    end
    return res
end

function heuristicSolve(
    graph;
    params=default_heur_params,
    timelimit::Float64=time() + 10.0,
)::StdResultWrapper

    node_multiplier::Ref{Float64} = Ref{Float64}(params.m_start)

    heur_graph = MetaGraph(
        DiGraph();
        label_type=Int64,
        vertex_data_type=HeurVertex,
        edge_data_type=HeurEdge,
        weight_function=heurWeight,
        default_weight=Inf,
        graph_data=graph[],
    )

    function createNode(node::Int64)::Nothing
        heur_graph[node] = HeurVertex(
            graph[node].is_s,  # is_s
            graph[node].is_t,  # is_t
            graph[node].p,  # p
            graph[node].ph,  # ph
            graph[node].p,  # weight
        )
        return nothing
    end

    function createEdge((n1, n2)::Tuple{Int64,Int64})::Nothing
        heur_graph[n1, n2] = HeurEdge(
            graph[n1, n2].d,  # d
            graph[n1, n2].big_d,  # big_d
            graph[n1, n2].d,  # weight
            Ref{HeurVertex}(heur_graph[n1]),  # src_node
            node_multiplier,
        )
        return nothing
    end

    # Prends un path au format MetaGraphsNext.jl
    function trueNodeWeights(path::Vector{Int64})::Dict{Int64,Float64}
        sorted_nodes::Vector{Int64} = sort(path, by=v -> heur_graph[v].ph, rev=true)
        res::Dict{Int64,Float64} = Dict{Int64,Float64}()
        d2::Float64 = heur_graph[].d2
        rem::Float64 = 0
        for node in sorted_nodes
            rem = min(d2, 2)
            res[node] = heur_graph[node].p + rem * heur_graph[node].ph
            d2 -= rem
        end
        return res
    end

    # Prends un path au format MetaGraphsNext.jl
    function trueEdgeWeights(path::Vector{Int64})::Dict{Tuple{Int64,Int64},Float64}
        arcs::Vector{Tuple{Int64,Int64}} = collect((path[i], path[i+1]) for i in 1:length(path)-1)
        sort!(arcs, by=((v1, v2),) -> heur_graph[v1, v2].d, rev=true)
        res::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()
        d1::Float64 = heur_graph[].d1
        rem::Float64 = 0
        for (v1, v2) in arcs
            rem = min(d1, heur_graph[v1, v2].big_d)
            res[(v1, v2)] = heur_graph[v1, v2].d * (1 + rem)
            d1 -= rem
        end
        return res
    end

    function updateNWeights(n_weights::Dict{Int64,Float64})::Nothing
        for (node, weight) in pairs(n_weights)
            heur_graph[node].weight = (
                params.v_smoother(n) * heur_graph[node].weight
                +
                (1 - params.v_smoother(n)) * weight
            )
        end
        return nothing
    end

    function updateAWeights(a_weights::Dict{Tuple{Int64,Int64},Float64})::Nothing
        for ((v1, v2), weight) in pairs(a_weights)
            heur_graph[v1, v2].weight = (
                params.e_smoother(n) * heur_graph[v1, v2].weight
                +
                (1 - params.e_smoother(n)) * weight
            )
        end
        return nothing
    end

    foreach(createNode, labels(graph))
    foreach(createEdge, edge_labels(graph))

    n::Int64 = 1
    best_sol_value::Float64 = Inf
    best_sol::Vector{Int64} = Vector{Int64}()
    while time() <= timelimit - 1
        state::Graphs.DijkstraState{Float64,Int64} = dijkstra_shortest_paths(
            heur_graph,
            code_for(heur_graph, heur_graph[].s),
        )
        path::Vector{Int64} = map(
            node -> label_for(heur_graph, node),  # Passe du format Graphs.jl au format MetaGraphsNext.jl
            mkPath(
                state;
                s=code_for(heur_graph, heur_graph[].s),
                t=code_for(heur_graph, heur_graph[].t),
            )
        )
        n_weights::Dict{Int64,Float64} = trueNodeWeights(path)
        updateNWeights(n_weights)

        a_weights::Dict{Tuple{Int64,Int64},Float64} = trueEdgeWeights(path)
        updateAWeights(a_weights)

        if sum(values(n_weights)) <= heur_graph[].big_s

            node_multiplier[] *= params.m_smoother(n)

            if sum(values(a_weights)) < best_sol_value
                best_sol_value = sum(values(a_weights))
                best_sol = path
            end

        else
            node_multiplier[] /= params.m_smoother(n)
        end
        n += 1
    end
    println(n)
    if best_sol_value < Inf
        return (
            is_feasible=true,
            proven_optimality=false,
            value=best_sol_value,
            lower_bound=0,
            upper_bound=best_sol_value,
            solution=best_sol,
        )
    else
        return (
            is_feasible=false,
            proven_optimality=false,
            value=typemax(Float64),
            lower_bound=0,
            upper_bound=typemax(Float64),
            solution=Vector{Int64}(),
        )
    end
end

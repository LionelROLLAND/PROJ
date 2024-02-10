using JuMP
using Gurobi
using Graphs
using MetaGraphsNext


function branch_and_cut(
    graph::MetaGraph;
    timelimit::Float64=time() + 3600.0,
)::ModResultWrapper
    s::Int64 = graph[].s
    t::Int64 = graph[].t

    eps::Float64 = 0.000_1
    function flow_creation(i::Int64)::Int64
        return Int64(graph[i].is_s) - Int64(graph[i].is_t)
    end

    modele_maitre = Model(Gurobi.Optimizer)
    set_silent(modele_maitre)
    @variable(modele_maitre, z >= 0)
    @variable(modele_maitre, a[edge_labels(graph)], Bin)  # vaut 1 ssi arc (i,j) choisi
    @objective(modele_maitre, Min, z)
    @constraint(modele_maitre, sum(graph[i, j].d * a[(i, j)] for (i, j) in edge_labels(graph)) <= z)
    # Conservation du flot :
    @constraint(
        modele_maitre,
        flow_cons[i in labels(graph)],
        (
            sum(a[(i, j)] for j in outneighbor_labels(graph, i))
            -
            sum(a[(j, i)] for j in inneighbor_labels(graph, i))
            ==
            flow_creation(i)
        ),
    )
    @constraint(
        modele_maitre,
        (
            sum(graph[i].p * a[(i, j)] for (i, j) in edge_labels(graph))
            +
            graph[graph[].t].p <= graph[].big_s
        ),
    )

    function heuristic_SP1(a_val)
        SP1_value = 0
        worst_weight::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()
        delta_d_heur::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()
        for (i, j) in edge_labels(graph)
            if (a_val[(i, j)] > 0)
                worst_weight[(i, j)] = a_val[(i, j)] * graph[i, j].d
                SP1_value += a_val[(i, j)] * graph[i, j].d
            end
            delta_d_heur[(i, j)] = 0
        end
        tuple_vec = sort(collect(worst_weight), by=x -> x.second, rev=true)

        d1::Float64 = graph[].d1
        rem::Float64 = 0
        for vec in tuple_vec
            i, j = vec.first
            rem = min(d1, graph[i, j].big_d)
            delta_d_heur[(i, j)] = rem
            SP1_value += a_val[(i, j)] * graph[i, j].d * rem
            d1 -= rem
            if d1 <= 0
                break
            end
        end
        return SP1_value, delta_d_heur
    end

    function heuristic_SP2(a_val)
        SP2_value = 0
        worst_weight::Dict{Int64,Float64} = Dict{Int64,Float64}()
        delta_p_heur::Dict{Int64,Float64} = Dict{Int64,Float64}()
        for (i, j) in edge_labels(graph)
            if (a_val[(i, j)] > 0)
                worst_weight[i] = graph[i].ph
                SP2_value += a_val[(i, j)] * graph[i].p
            end
            delta_p_heur[i] = 0
        end
        worst_weight[t] = graph[t].ph
        SP2_value += graph[t].p
        delta_p_heur[t] = 0
        tuple_vec = sort(collect(worst_weight), by=x -> x.second, rev=true)

        d2::Float64 = graph[].d2
        rem::Float64 = 0
        for vec in tuple_vec
            i = vec.first
            rem = min(d2, 2)
            delta_p_heur[i] = rem
            if i == t
                SPA_value += graph[i].ph * rem
            else
                SP2_value += sum(a_val[(i, j)] for j in outneighbor_labels(graph, i); init=0) * graph[i].ph * rem
            end
            d2 -= rem
            if d2 <= 0
                break
            end
        end
        return SP2_value, delta_p_heur
    end

    function my_lazy_callback(cb_data)
        # On récupère la valeur de a
        a_val = Dict((i, j) => callback_value(cb_data, a[(i, j)]) for (i, j) in edge_labels(graph))
        z_val = callback_value(cb_data, z)

        # Résolution SP1
        SP1_value, delta_d_sol = heuristic_SP1(a_val)
        if SP1_value > z_val + eps
            cstr = @build_constraint(
                sum(graph[i, j].d * (1 + delta_d_sol[(i, j)]) * a[(i, j)] for (i, j) in edge_labels(graph)) <= z)
            MOI.submit(modele_maitre, MOI.LazyConstraint(cb_data), cstr)
            println("Added SP1 constraint")
        end

        # Résolution SP2
        SP2_value, delta_p_sol = heuristic_SP2(a_val)
        if SP2_value > graph[].big_s + eps
            cstr = @build_constraint(
                sum((graph[i].p + delta_p_sol[i] * graph[i].ph) * a[(i, j)] for (i, j) in edge_labels(graph)) +
                graph[t].p + delta_p_sol[t] * graph[t].ph <= graph[].big_s
            )
            MOI.submit(modele_maitre, MOI.LazyConstraint(cb_data), cstr)
            println("Added SP2 constraint")
        end
    end

    start = time()
    MOI.set(modele_maitre, MOI.LazyConstraintCallback(), my_lazy_callback)
    set_attribute(modele_maitre, "TimeLimit", timelimit - time())
    JuMP.optimize!(modele_maitre)
    println("elapsed time : ", time() - start)
    prim_stat::ResultStatusCode = primal_status(modele_maitre)
    a_sol::Dict{Tuple{Int64,Int64},Float64} = (
        prim_stat == FEASIBLE_POINT ?
        Dict{Tuple{Int64,Int64},Float64}((i, j) => value(a[(i, j)]) for (i, j) in edge_labels(graph)) :
        Dict{Tuple{Int64,Int64},Float64}((i, j) => 0.0 for (i, j) in edge_labels(graph))
    )
    return (
        primal_status=prim_stat,
        dual_status=dual_status(modele_maitre),
        term_status=termination_status(modele_maitre),
        obj_value=(prim_stat == FEASIBLE_POINT ? objective_value(modele_maitre) : Inf),
        lower_bound=objective_bound(modele_maitre),
        upper_bound=(prim_stat == FEASIBLE_POINT ? objective_value(modele_maitre) : Inf),
        a=a_sol,
    )
end

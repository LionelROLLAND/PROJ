using Graphs
using MetaGraphsNext
using JuMP

ModResultWrapper = @NamedTuple begin
    primal_status::ResultStatusCode
    dual_status::ResultStatusCode
    term_status::TerminationStatusCode
    obj_value::Float64
    lower_bound::Float64
    upper_bound::Float64
    a::Dict{Tuple{Int64,Int64},Float64}
end

StdResultWrapper = @NamedTuple begin
    is_feasible::Bool
    proven_optimality::Bool
    value::Float64
    lower_bound::Float64
    upper_bound::Float64
    solution::Vector{Int64}
end

RawData = @NamedTuple begin
    n::Int64
    s::Int64
    t::Int64
    big_s::Int64
    d1::Int64
    d2::Int64
    p::Vector{Int64}
    ph::Vector{Int64}
    arcs::Dict{Tuple{Int64,Int64},NamedTuple{(:d, :big_d),Tuple{Int64,Float64}}}
end


function readInstance(fd::IOStream)::RawData
    lines::Vector{String} = readlines(fd)
    n::Int64 = parse(Int64, split(lines[1], " ")[end])
    s::Int64 = parse(Int64, split(lines[2], " ")[end])
    t::Int64 = parse(Int64, split(lines[3], " ")[end])
    big_s::Int64 = parse(Int64, split(lines[4], " ")[end])
    d1::Int64 = parse(Int64, split(lines[5], " ")[end])
    d2::Int64 = parse(Int64, split(lines[6], " ")[end])
    str_p::Vector{String} = split(
        strip(lines[7], Set{Char}(['p', '=', ' ', '[', ']', ';'])),
        ", ";
        keepempty=false,
    )
    p::Vector{Int64} = collect(parse(Int64, str_i) for str_i in str_p)
    str_ph::Vector{String} = split(
        strip(lines[8], Set{Char}(['p', 'h', '=', ' ', '[', ']', ';'])),
        ", ";
        keepempty=false,
    )
    ph::Vector{Int64} = collect(parse(Int64, str_i) for str_i in str_ph)
    arc_dict::Dict{Tuple{Int64,Int64},NamedTuple{(:d, :big_d),Tuple{Int64,Float64}}} = Dict()
    for line in lines[10:end]
        str_arc::Vector{String} = split(
            strip(line, Set{Char}([' ', ';', '[', ']'])),
            " ";
            keepempty=false,
        )
        i::Int64 = parse(Int64, str_arc[1])
        j::Int64 = parse(Int64, str_arc[2])
        d::Int64 = parse(Int64, str_arc[3])
        big_d::Float64 = parse(Float64, str_arc[4])
        arc_dict[(i, j)] = (d=d, big_d=big_d)
    end
    return (n=n, s=s, t=t, big_s=big_s, d1=d1, d2=d2, p=p, ph=ph, arcs=arc_dict)
end

function weight_fun(
    ntup::NamedTuple{
        (:d, :big_d, :weight),
        Tuple{Int64,Float64,Ref{Float64}}
    }
)::Float64
    return ntup.weight[]
end

function graphFromData(data::RawData)::MetaGraph
    if data.s == data.t
        error("The source and the sink are the same node.")
    end
    graph = MetaGraph(
        DiGraph();
        label_type=Int64,
        vertex_data_type=NamedTuple{(:is_s, :is_t, :p, :ph),Tuple{Bool,Bool,Int64,Int64}},
        edge_data_type=NamedTuple{(:d, :big_d, :weight),Tuple{Int64,Float64,Ref{Float64}}},
        weight_function=weight_fun,
        default_weight=Inf,
        graph_data=(s=data.s, t=data.t, big_s=data.big_s, d1=data.d1, d2=data.d2),
    )
    foreach(i -> graph[i] = (is_s=(i == data.s), is_t=(i == data.t), p=data.p[i], ph=data.ph[i]), 1:data.n)
    for ((i, j), arc) in data.arcs
        # On degage les arcs qui entrent en s ou sortent de t, on en a pas besoin et ca complexifie les
        # formulations pour rien
        if !(graph[i].is_t) && !(graph[j].is_s)
            graph[i, j] = (d=arc.d, big_d=arc.big_d, weight=Ref{Float64}(Inf))
        end
    end
    return graph
end

function mkPath(a::Dict{Tuple{Int64,Int64},Float64})::Vector{Int64}
    building_dict::Dict{Int64,Int64} = Dict{Int64,Int64}()
    eps = 0.000_1
    for ((i, j), val) in pairs(a)
        int_val::Int64 = round(Int64, val)
        abs(val - round(val)) <= eps && int_val in [0, 1] ? nothing : error("a is not integer feasible.")
        if int_val == 1
            haskey(building_dict, i) ? error("Two arcs coming from the same node.") : building_dict[i] = j
        end
    end
    current_i::Int64 = only(setdiff(keys(building_dict), values(building_dict)))
    path::Vector{Int64} = [current_i]
    while haskey(building_dict, current_i)
        current_i = building_dict[current_i]
        push!(path, current_i)
    end
    return path
end

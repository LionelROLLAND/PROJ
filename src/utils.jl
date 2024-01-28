using Graphs
using MetaGraphsNext


function readInstance(fd::IOStream)::NamedTuple{
    (:n, :s, :t, :big_s, :d1, :d2, :p, :ph, :arcs),
    Tuple{
        Int64,
        Int64,
        Int64,
        Int64,
        Int64,
        Int64,
        Vector{Int64},
        Vector{Int64},
        Dict{Tuple{Int64,Int64},NamedTuple{(:d, :big_d),Tuple{Int64,Float64}}},
    },
}
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

function graphFromData(
    data::NamedTuple{
        (:n, :s, :t, :big_s, :d1, :d2, :p, :ph, :arcs),
        Tuple{
            Int64,
            Int64,
            Int64,
            Int64,
            Int64,
            Int64,
            Vector{Int64},
            Vector{Int64},
            Dict{Tuple{Int64,Int64},NamedTuple{(:d, :big_d),Tuple{Int64,Float64}}},
        },
    }
)::MetaGraph
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
    for i in 1:data.n
        if i == data.s
            graph[i] = (is_s=true, is_t=false, p=data.p[i], ph=data.ph[i])
        elseif i == data.t
            graph[i] = (is_s=false, is_t=true, p=data.p[i], ph=data.ph[i])
        else
            graph[i] = (is_s=false, is_t=false, p=data.p[i], ph=data.ph[i])
        end
    end
    for ((i, j), arc) in data.arcs
        graph[i, j] = (d=arc.d, big_d=arc.big_d, weight=Ref{Float64}(Inf))
    end
    return graph
end

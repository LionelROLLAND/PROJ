using Plots

include("latex.jl")

function processResultsForPerfDiag(
    src_files::Vector{String};
    not_robust_method::String,
)::Dict{String,Vector{Float64}}
    eps::Float64 = 0.001
    preprocessed::ProcessedResultsForTable = processResultsForTable(
        src_files;
        not_robust_method=not_robust_method,
    )
    s_time_lists::Dict{String,Vector{Float64}}
    robust_list::Vector{InstanceResEntry} = preprocessed.robust_list
    lower_bound::Dict{String,Float64} = preprocessed.lower_bounds
    for (instance, res) in robust_list
        for (method, sol) in pairs(res)
            if sol.is_feasible && sol.value <= lower_bound[instance] + eps
                s_times::Vector{Float64} = get!(s_time_lists, method, [])
                push!(s_times, sol.solving_time)
            end
        end
    end
    foreach(s_time -> sort!(s_time), values(s_time_lists))
    return s_time_lists
end

function savePerfDiag(
    s_time_lists::Dict{String,Vector{Float64}};
    until::Float64,
    size::Tuple{Int64,Int64}=(1000, 1000),
)::Nothing
    points::Dict{String,Vector{Tuple{Float64,Int64}}} = Dict{String,Vector{Tuple{Float64,Int64}}}(
        method => [(0.0, 0)] for method in keys(s_time_lists)
    )
    for (method, s_times) in pairs(s_time_lists)
        last_i::Int64 = 0
        for (i, val) in enumerate(s_time_lists[method])
            if val > until
                break
            end
            push!(s_times, (val, i - 1), (val, i))
            last_i = i
        end
        push!(s_times, (until, last_i))
    end
    gr(size=size) # Not finished
end

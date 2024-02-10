using Graphs
using MetaGraphsNext
import JSON

include("utils.jl")

function testMethod(;
    method::Function,
    method_name::String,
    save::String,
    instance_dir::String,
)::Nothing
    save_file = ReentrantLock()
    function instanceWrapper(filename::String)::Nothing
        println("Processing $filename...")
        graph::MetaGraph = open(
            graphFromData ∘ readInstance,
            joinpath(instance_dir, filename),
        )
        elapsed = @elapsed result::StdResultWrapper = method(graph)
        result_dict::Dict{String,Union{Int64,Float64,String,Vector{Int64}}} = Dict{
            String,Union{Int64,Float64,String,Vector{Int64}}
        }(
            "file" => filename,
            "method" => method_name,
            "is_feasible" => Int64(result.is_feasible),
            "proven_optimality" => Int64(result.proven_optimality),
            "value" => result.value,
            "lower_bound" => result.lower_bound,
            "upper_bound" => result.upper_bound,
            "solution" => result.solution,
            "solving_time" => elapsed,
        )
        while !trylock(save_file)
            sleep(0.01)
        end
        data = open(fd -> JSON.parse(fd), save, read=true)
        push!(data, result_dict)
        open(fd -> JSON.print(fd, data), save, write=true)
        unlock(save_file)
        return nothing
    end

    already_done::Set{String} = Set{String}()
    if !isfile(save)
        open(
            fd -> JSON.print(fd, []),
            save,
            write=true,
            create=true,
        )
    else
        union!(
            already_done,
            Iterators.map(
                entry -> entry["file"],
                Iterators.filter(
                    entry -> entry["method"] == method_name,
                    open(fd -> JSON.parse(fd), save, read=true),
                ),
            )
        )
    end
    @sync begin
        for file in Iterators.filter(f -> !(f in already_done), readdir(instance_dir))
            Threads.@spawn instanceWrapper($file)
        end
    end
    return nothing
end

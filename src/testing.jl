using Graphs
using MetaGraphsNext
import JSON

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
        elapsed = @elapsed result = method(graph)
        result_dict::Dict{String,Union{String,Int64,Float64}} = Dict{String,Union{String,Int64,Float64}}(
            "file" => filename,
            "method" => method_name,
            "is_feasible" => Int64(result.is_feasible),
            "proven_optimality" => Int64(result.proven_optimality),
            "value" => result.value,
            "lower_bound" => result.lower_bound,
            "upper_bound" => result.upper_bound,
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

    if !isfile(save)
        open(
            fd -> JSON.print(fd, []),
            save,
            write=true,
            create=true,
        )
    end
    save_file = ReentrantLock()
    @sync begin
        for file in readdir(instance_dir)
            Threads.@spawn instanceWrapper(file)
        end
    end
    return nothing
end

MethodResult = @NamedTuple begin
    is_feasible::Bool
    proven_optimality::Bool
    value::Float64
    lower_bound::Float64
    upper_bound::Float64
    solving_time::Float64
end

InstanceResEntry = @NamedTuple begin
    instance::String
    res::Dict{String,MethodResult}
end

ProcessedResultsForTable = @NamedTuple begin
    robust_list::Vector{InstanceResEntry}
    not_robust_data::Dict{String,Dict{String,MethodResult}}
    lower_bounds::Dict{String,Float64}
    upper_bounds::Dict{String,Float64}
end

function nameSorter(s::String)::Tuple{Int64,String}
    n_str, suff = split(s, "_")
    return parse(Int64, n_str), suff
end

function gap(val::Float64, ref::Float64)::Float64
    if ref == 0
        if val == 0
            return 0
        else
            return Inf
        end
    else
        return (val - ref) / ref
    end
end

function percentify(val::Float64)::String
    if val == Inf
        return raw"\infty"
    else
        return string(round(100 * val, digits=1)) * raw" \ \% "
    end
end

function minutify(seconds::Float64)::String
    mins, secs = divrem(seconds, 60)
    if mins == 0
        return "$(round(secs, RoundUp; digits=1))s"
    else
        return "$(Int64(mins))m$(round(secs, RoundUp; digits=1))s"
    end
end

function latexify(s::String)::String
    return replace(s, "_" => raw"\_")
end


function processResultsForTable(
    src_files::Vector{String};
    not_robust_method::String,
)::ProcessedResultsForTable
    # Dict(instance => Dict(method => result))
    robust_data::Dict{String,Dict{String,MethodResult}} = Dict{String,Dict{String,MethodResult}}()
    not_robust_data::Dict{String,Dict{String,MethodResult}} = Dict{String,Dict{String,MethodResult}}()
    for src_file in src_files
        src_data = open(fd -> JSON.parse(fd), src_file, read=true)
        for entry_dict in src_data
            instance_entry::Dict{String,MethodResult} = get!(
                entry_dict["method"] == not_robust_method ? not_robust_data : robust_data,
                entry_dict["file"],
                Dict{String,MethodResult}(),
            )
            if haskey(instance_entry, entry_dict["method"])
                @warn begin
                    "Several resolutions found with the same method ($(entry_dict["method"]))" *
                    " for the same instance ($(entry_dict["file"]))."
                end
            end

            instance_entry[entry_dict["method"]] = (
                is_feasible=Bool(entry_dict["is_feasible"]),
                proven_optimality=Bool(entry_dict["proven_optimality"]),
                value=entry_dict["value"],
                lower_bound=entry_dict["lower_bound"],
                upper_bound=entry_dict["upper_bound"],
                solving_time=entry_dict["solving_time"],
            )

        end
    end

    lower_bound::Dict{String,Float64} = Dict{String,Float64}()
    upper_bound::Dict{String,Float64} = Dict{String,Float64}()
    for (instance, instance_entry) in pairs(robust_data)
        lower_bound[instance] = maximum(
            method -> method.lower_bound,
            values(instance_entry),
            init=0,
        )
        upper_bound[instance] = minimum(
            method -> method.upper_bound,
            values(instance_entry),
            init=Inf,
        )
    end
    instance_list::Vector{InstanceResEntry} = collect(
        (instance=instance, res=entry) for (instance, entry) in pairs(robust_data)
    )
    sort!(instance_list; by=ins_res_entry -> nameSorter(ins_res_entry.instance))
    return (
        robust_list=instance_list,
        not_robust_data=not_robust_data,
        lower_bounds=lower_bound,
        upper_bounds=upper_bound,
    )
end


function writeTable(
    methods::Vector{String};
    results::ProcessedResultsForTable,
    template_file::String,
    output::IOStream,
)::Nothing

    robust_list = results.robust_list
    not_robust_data = results.not_robust_data
    lower_bound = results.lower_bounds
    upper_bound = results.upper_bounds

    function insertTable(io::IOStream)::Nothing
        println(io, raw"\begin{tabular}{", join(fill("c", 2 + 2 * length(methods)), " "), "}")
        println(io, raw"\hline")
        println(
            io,
            raw"\multirow{2}{*}{Instance} & \multirow{2}{*}{PR} & ",
            join(collect("\\multicolumn{2}{c }{$method}" for method in methods), " & "),
            "\\\\",
        )
        println(
            io,
            "& & ",
            join(fill("Time & Gap", length(methods)), " & "),
            "\\\\",
        )
        println(io, raw"\hline")
        eps = 0.001
        for (instance, entry) in robust_list
            print(io, latexify(splitext(instance)[1]), " & ")

            lower_exact_value::Float64 = 0
            if haskey(not_robust_data, instance)
                meth_entry = not_robust_data[instance] |> values |> only
                lower_exact_value = meth_entry.lower_bound
                if lower_bound[instance] < lower_exact_value
                    lower_bound[instance] = lower_exact_value
                end
            end
            upper_exact_value::Float64 = minimum(
                meth -> meth.upper_bound,
                (
                    values(entry),
                    get(not_robust_data, instance, Dict{String,MethodResult}()) |> values,
                ) |> Iterators.flatten;
                init=Inf,
            )
            if !haskey(not_robust_data, instance)
                print(io, "? ")
            elseif upper_bound[instance] * upper_exact_value <= lower_exact_value * lower_bound[instance] + eps
                print(
                    io,
                    raw"$",
                    gap(upper_bound[instance], lower_exact_value) |> percentify,
                    raw"$ ",
                )
            else
                print(
                    io,
                    raw"$",
                    gap(lower_bound[instance], upper_exact_value) |> percentify,
                    raw" \leq \quad \leq ",
                    gap(upper_bound[instance], lower_exact_value) |> percentify,
                    raw"$ ",
                )
            end
            for method in methods
                print(io, "& ")
                if haskey(entry, method)
                    meth_entry = entry[method]
                    print(io, minutify(meth_entry.solving_time), " & ")
                    if meth_entry.is_feasible
                        print(io, raw"$", gap(meth_entry.value, lower_bound[instance]) |> percentify, raw"$ ")
                    else
                        print(io, "- ")
                    end
                else
                    print(io, " & ")
                end
            end
            println(io, "\\\\")
        end
        println(io, raw"\hline")
        println(io, raw"\end{tabular}")
        return nothing
    end

    pref, suff = open(
        fd -> split(read(fd, String), "%insert table%"),
        template_file,
        read=true,
    )
    print(output, pref)
    insertTable(output)
    print(output, suff)
    return nothing
end

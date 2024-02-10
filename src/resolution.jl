using JuMP
using Gurobi
using Graphs
using MetaGraphsNext


function resolution(
    graph::MetaGraph;
    timelimit::Float64=time() + 3600.0,
)::ModResultWrapper

    function flow_creation(i::Int64)::Int64
        return Int64(graph[i].is_s) - Int64(graph[i].is_t)
    end

    m = Model(Gurobi.Optimizer)
    set_silent(m)
    @variable(m, a[edge_labels(graph)], Bin)  # vaut 1 ssi arc (i,j) choisi
    @objective(m, Min, sum(graph[i, j].d * a[(i, j)] for (i, j) in edge_labels(graph)))
    # Conservation du flot :
    @constraint(
        m,
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
        m,
        (
            sum(graph[i].p * a[(i, j)] for (i, j) in edge_labels(graph))
            +
            graph[graph[].t].p <= graph[].big_s
        ),
    )

    start = time()
    s::Int64 = graph[].s
    t::Int64 = graph[].t
    set_attribute(m, "TimeLimit", timelimit - time())
    JuMP.optimize!(m)
    m_value = JuMP.objective_value(m)
    println("Objective value m: ", m_value)
    println("elapsed time : ", time() - start)
    return (
        primal_status=primal_status(m),
        dual_status=dual_status(m),
        term_status=termination_status(m),
        obj_value=objective_value(m),
        lower_bound=objective_bound(m),
        upper_bound=(primal_status(m) == FEASIBLE_POINT ? objective_value(m) : Inf),
        a=Dict{Tuple{Int64,Int64},Float64}((i, j) => value(a[(i, j)]) for (i, j) in edge_labels(graph)),
    )
end

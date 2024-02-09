using JuMP
using Gurobi

include("utils.jl")

function dualSolve(
    graph::MetaGraph;
    timelimit::Float64=time() + 3600.0,
)::ModResultWrapper
    s::Int64 = graph[].s
    t::Int64 = graph[].t
    m = Model(Gurobi.Optimizer)
    set_silent(m)
    function flow_creation(i::Int64)::Int64
        return Int64(graph[i].is_s) - Int64(graph[i].is_t)
    end
    @variable(m, a[edge_labels(graph)], Bin)
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
        )
    )
    @variable(m, lambda_d >= 0)
    @variable(m, mu_d[edge_labels(graph)] >= 0)
    # Contrainte de (D_d)
    @constraint(
        m,
        d_dual_cons[(i, j) in edge_labels(graph)],
        lambda_d + mu_d[(i, j)] >= graph[i, j].d * a[(i, j)],
    )
    @variable(m, lambda_p >= 0)
    @variable(m, mu_p[labels(graph)] >= 0)
    # Contrainte de (D_p)
    @constraint(
        m,
        p_dual_cons[i in labels(graph)],
        lambda_p + mu_p[i]
        >=
        graph[i].ph * (graph[i].is_t + sum(a[(i, j)] for j in outneighbor_labels(graph, i))),
    )
    # Contrainte correspondant a l'objectif de (D_p)
    @constraint(
        m,
        p_dual_obj,
        (
            graph[t].p
            +
            sum(graph[i].p * a[(i, j)] for i in labels(graph) for j in outneighbor_labels(graph, i))
            +
            2 * sum(mu_p) + lambda_p * graph[].d2
            <=
            graph[].big_s
        )
    )
    # Objectif
    @objective(
        m,
        Min,
        graph[].d1 * lambda_d
        +
        sum(
            graph[i, j].d * a[(i, j)]
            +
            graph[i, j].big_d * mu_d[(i, j)]
            for (i, j) in edge_labels(graph)
        )
    )
    set_attribute(m, "TimeLimit", timelimit - time())
    optimize!(m)
    ub::Float64 = Inf
    prim_stat::ResultStatusCode = primal_status(m)
    if prim_stat == FEASIBLE_POINT
        ub = objective_value(m)
    end
    return (
        primal_status=primal_status(m),
        dual_status=dual_status(m),
        term_status=termination_status(m),
        obj_value=objective_value(m),
        lower_bound=objective_bound(m),
        upper_bound=ub,
        a=Dict{Tuple{Int64,Int64},Float64}((i, j) => value(a[(i, j)]) for (i, j) in edge_labels(graph)),
    )
end

function completeModelWrapper(method::Function; time_budget::Float64=60.0)::Function
    function wrapped(graph::MetaGraph)::StdResultWrapper
        result::ModResultWrapper = method(graph; timelimit=time() + time_budget)
        return (
            is_feasible=(result.primal_status == FEASIBLE_POINT),
            proven_optimality=(result.term_status == OPTIMAL),
            value=result.obj_value,
            lower_bound=result.lower_bound,
            upper_bound=result.upper_bound,
            solution=(result.primal_status == FEASIBLE_POINT ? mkPath(result.a) : []),
        )
    end
    return wrapped
end

using Graphs
using MetaGraphsNext
using JuMP
using Gurobi

function plans_coupants(
    graph::MetaGraph;
    timelimit::Float64=time() + 3600.0,
)::ModResultWrapper

    s::Int64 = graph[].s
    t::Int64 = graph[].t

    function flow_creation(i::Int64)::Int64
        return Int64(graph[i].is_s) - Int64(graph[i].is_t)
    end

    function heuristic_SP1(a_val)
        SP1_value = 0
        worst_weight::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()
        delta_d_heur::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()
        for (i, j) in edge_labels(graph)
            if (a_val[(i, j)] > 0)
                worst_weight[(i, j)] = graph[i, j].d * graph[i, j].big_d
                SP1_value += graph[i, j].d
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
            SP1_value += graph[i, j].d * rem
            d1 -= rem
            if (d1 == 0)
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
                worst_weight[i] = 2 * graph[i].ph
                SP2_value += graph[i].p
            end
            delta_p_heur[i] = 0
        end
        worst_weight[t] = 2 * graph[t].ph
        SP2_value += graph[t].p
        delta_p_heur[t] = 0
        tuple_vec = sort(collect(worst_weight), by=x -> x.second, rev=true)

        d2::Float64 = graph[].d2
        rem::Float64 = 0
        for vec in tuple_vec
            i = vec.first
            rem = min(d2, 2)
            delta_p_heur[i] = rem
            SP2_value += graph[i].ph * rem
            d2 -= rem
            if (d2 == 0)
                break
            end
        end
        return SP2_value, delta_p_heur
    end

    start = time()
    # __________MODELE_MAITRE__________
    #declaration du modele
    modele_maitre = Model(Gurobi.Optimizer)
    modele_maitre = Model(
        optimizer_with_attributes(
            Gurobi.Optimizer, "Presolve" => 0
        )
    )
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
    cpt = 0
    best_sup::Float64 = Inf
    SP1_value::Float64 = 0.0
    SP2_value::Float64 = 0.0
    Z::Float64 = 0
    while time() < timelimit
        JuMP.optimize!(modele_maitre)
        m_value = JuMP.objective_value(modele_maitre)
        Z = JuMP.value(z)
        A = JuMP.value.(a)
        for (i, j) in edge_labels(graph)
            if A[(i, j)] > 0
                println("(", i, ",", j, ") = ", A[(i, j)])
            end
        end
        println("Objective value modele_maitre: ", m_value)
        cpt += 1
        # # Résolution SP1
        SP1_value, delta_d_sol = heuristic_SP1(A)
        # # Résolution SP2
        SP2_value, delta_p_sol = heuristic_SP2(A)
        if (round(SP1_value - Z; digits=4) <= 0 && SP2_value <= graph[].big_s)
            break
        end
        println("Objective value SP1: ", SP1_value)
        if (round(SP1_value - Z; digits=4) > 0)
            @constraint(modele_maitre,
                sum(graph[i, j].d * (1 + delta_d_sol[(i, j)]) * a[(i, j)] for (i, j) in edge_labels(graph)) <= z)
            if (SP1_value < best_sup)
                best_sup = SP1_value
            end
            @constraint(modele_maitre, z <= best_sup)
            println("Violée 1")
        end

        println("Objective value SP2: ", SP2_value)
        if (SP2_value > graph[].big_s)
            @constraint(modele_maitre,
                sum((graph[i].p + delta_p_sol[i] * graph[i].ph) * a[(i, j)] for (i, j) in edge_labels(graph)) +
                graph[t].p + delta_p_sol[t] * graph[t].ph <= graph[].big_s)
            println("Violée 2")
        end
    end
    is_feasible::Bool = false
    is_optimal::Bool = false
    if round(SP2_value - graph[].big_s; digits=4) <= 0
        is_feasible = true
        if round(SP1_value - Z; digits=4) <= 0
            is_optimal = true
        end
    end
    println("COMPTEUR = ", cpt)
    println("elapsed time : ", time() - start)
    return (
        primal_status=(is_feasible ? FEASIBLE_POINT : INFEASIBLE_POINT),
        dual_status=dual_status(modele_maitre),
        term_status=(is_optimal ? OPTIMAL : TIME_LIMIT),
        obj_value=objective_value(modele_maitre),
        lower_bound=objective_bound(modele_maitre),
        upper_bound=(is_feasible ? objective_value(modele_maitre) : Inf),
        a=Dict{Tuple{Int64,Int64},Float64}((i, j) => value(a[(i, j)]) for (i, j) in edge_labels(graph)),
    )
end

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
                SP2_value += graph[i].ph * rem
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

    start = time()
    # __________MODELE_MAITRE__________
    #declaration du modele
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
    cpt = 0

    best_prim_stat::ResultStatusCode = NO_SOLUTION
    best_dual_stat::ResultStatusCode = NO_SOLUTION
    best_term_stat::TerminationStatusCode = OPTIMIZE_NOT_CALLED
    best_val::Float64 = Inf
    best_ub::Float64 = Inf
    best_lb::Float64 = 0
    best_a::Dict{Tuple{Int64,Int64},Float64} = Dict{Tuple{Int64,Int64},Float64}()

    eps = 0.000_1
    while time() < timelimit && best_term_stat != OPTIMAL
        set_attribute(modele_maitre, "TimeLimit", timelimit - time())
        JuMP.optimize!(modele_maitre)
        if best_term_stat == OPTIMIZE_NOT_CALLED
            best_term_stat = TIME_LIMIT
        end
        if primal_status(modele_maitre) != FEASIBLE_POINT
            break
        end
        best_lb = m_value = JuMP.objective_value(modele_maitre)
        Z = JuMP.value(z)
        A = Dict{Tuple{Int64,Int64},Float64}(
            (i, j) => JuMP.value(a[(i, j)])
            for (i, j) in edge_labels(graph)
        )
        println("Objective value modele_maitre: ", m_value)
        cpt += 1
        # # Résolution SP1
        SP1_value, delta_d_sol = heuristic_SP1(A)
        # # Résolution SP2
        SP2_value, delta_p_sol = heuristic_SP2(A)

        println("Objective value SP2: ", SP2_value)
        if SP2_value > graph[].big_s + eps
            @constraint(
                modele_maitre,
                (
                    sum(
                        (graph[i].p + delta_p_sol[i] * graph[i].ph) * a[(i, j)]
                        for (i, j) in edge_labels(graph)
                    )
                    +
                    graph[t].p + delta_p_sol[t] * graph[t].ph
                    <=
                    graph[].big_s
                ),
            )
            println("Violée 2")
        else
            # La solution est realisable
            if SP1_value <= best_ub
                # Le <= (et pas <) est hyper important ici a cause de la mise a jour de
                # best_term_stat
                best_prim_stat = FEASIBLE_POINT
                best_dual_stat = UNKNOWN_RESULT_STATUS
                best_val = SP1_value
                if SP1_value < best_ub
                    best_ub = SP1_value
                    @constraint(modele_maitre, z <= best_ub)
                end
                best_a = deepcopy(A)
                if SP1_value <= Z + eps
                    best_term_stat = OPTIMAL
                end
            end
        end

        println("Objective value SP1: ", SP1_value)
        if SP1_value > Z + eps
            @constraint(
                modele_maitre,
                sum(
                    graph[i, j].d * (1 + delta_d_sol[(i, j)]) * a[(i, j)]
                    for (i, j) in edge_labels(graph)
                ) <= z
            )
            println("Violée 1")
        end

    end
    println("COMPTEUR = ", cpt)
    println("elapsed time : ", time() - start)
    return (
        primal_status=best_prim_stat,
        dual_status=best_dual_stat,
        term_status=best_term_stat,
        obj_value=best_val,
        lower_bound=best_lb,
        upper_bound=best_ub,
        a=best_a,
    )
end

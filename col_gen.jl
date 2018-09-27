# This file matches all data and creates the relevant dataframe
using CSV
using DataFrames
using CPLEX
using JuMP
using LightGraphs
using SimpleWeightedGraphs
using MetaGraphs
using GraphLayout
using Combinatorics
using OffsetArrays

include("helper.jl")

struct soln
    edge
    vertex
end

if length(ARGS) >= 1
    println(ARGS[1])
    global n = parse(Int64,ARGS[1])
elseif !isdefined(:n)
    global n = 10
end

# max_weight = Int(round(n/2))
max_weight = 3

positions = CSV.read("./data/positions_"*string(n)*".csv")
positions[:weight] = ones(n)
positions[:profit] = positions[:profit]*3

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')./2
end

distances = get_distance(positions)

# if n ==5
#     distances = [0.0 3 5 7 8;
#                  3 0 6 6 3;
#                  5 6 0 3 5;
#                  7 6 3 0 2;
#                  8 3 5 2 0]
#     max_weight = 2
#     positions[:profit] = [0,12,8,8,3]
#     positions[:weight] = [0,1,1,1,1]
# end
################################
# Create the RMP_τ^(I^P)
################################
global τ = 1
println("Create the RMP")
ϵ_α = 100000
ϵ_β = 100000
ϵ_μ = 10

α_tilde = fill!(OffsetArray{Float64}(3:n),0)
β_tilde = 0

I_P = Set{Int64}()
v_lb = -Inf
v_ub = Inf
c_1 = Inf

columns = Dict()

rmp_τ = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 0, CPXPARAM_Preprocessing_Dual=-1, CPXPARAM_MIP_Display = 3))

@variable(rmp_τ, η_ub[3:n] >=0)
@variable(rmp_τ, η_lb[3:n] >=0)
@variable(rmp_τ, η_lb_knap >=0)
@variable(rmp_τ, η_ub_knap >=0)

# Add the objective function
@objective(rmp_τ, Min, sum( (α_tilde[v] + ϵ_α) * η_ub[v] - (α_tilde[v] - ϵ_α) * η_lb[v] for v=3:n) +
    (β_tilde + ϵ_β)*η_ub_knap - (β_tilde - ϵ_β)*η_lb_knap)

# Add constraint (10a)
@constraint(rmp_τ, vertex_constraints[i=3:n], η_ub[i] - η_lb[i] == 0 )
@constraint(rmp_τ, knapsack_constraint, η_ub_knap - η_lb_knap <= 0 )

# Create the constraint refs
constraint_refs = Vector{ConstraintRef}(0)
for i=3:n
    push!(constraint_refs, vertex_constraints[i])
end
push!(constraint_refs, knapsack_constraint)


function create_reduced_costs(rmp)
    α = getdual(vertex_constraints)
    β = getdual(knapsack_constraint)
    w_hat = copy(convert(Array{Float64,2},distances))
    for i=3:n
        for j=i+1:n
            w_hat[i,j] -= α[i] + α[j]
            w_hat[j,i] = w_hat[i,j]
        end
    end

    r = positions[:profit]
    r_hat = copy(convert(Array{Float64,1},r))
    for v=3:n
        r_hat[v] = r[v] - 2*α[v] + positions[:weight][v] * β
    end

    return w_hat, r_hat
end

function create_cgp(rmp)
    w_hat, r_hat = create_reduced_costs(rmp)
    cgp = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))

    @variable(cgp, edge[i=1:n, j=i+1:n], Bin)
    @variable(cgp, vertex[i=1:n], Bin)
    @objective(cgp, Min, sum(w_hat[i,j]*edge[i,j] for i=1:n, j=i+1:n) - sum(r_hat[i] * vertex[i] for i=2:n))

    # Constraint (6c)
    @constraint(cgp, sum(edge[1,i] for i=2:n) == 2)
    # Constraint (2a)
    @constraint(cgp, sum(edge[i,j] for i=1:n, j=i+1:n) == sum(vertex[i] for i=1:n))
    @constraint(cgp, vertex[1] == 1)
    # Constraints (2b)
    println("Creating the generalised subtour elimination constraints")
    # TODO: Fix this constraint, edges must have one vertex IN and one vertex OUT of the set
    for num in 1:n
        println("Up to $num cities")
        for list in combinations(2:n, num)
            println(list)
            println(length(list))
            @constraint(cgp, sum(edge[i,j] for i in list, j in list if i < j) >= 1 + sum(vertex[i] for i in list if i!=1) - length(list))
        end
    end
    return cgp
end

function swap_edges_remove!(vertex, tree, edges)
    # First find the best two edges to swap
    best_swap_value = Inf
    best_swap = [(1,2),(2,3)]
    function find_weight(edge, edge_df)
        if edge[1] < edge[2]
            return edge_df[(edge_df[:in].==edge[1]) .& (edge_df[:out].==edge[2]),:weight][1]
        elseif edge[1] > edge[2]
            return edge_df[(edge_df[:in].==edge[2]) .& (edge_df[:out].==edge[1]),:weight][1]
        else
            println("ERROR")
        end
    end
    incident_edge_list = outneighbors(tree,vertex)
    # println(incident_edge_list)
    for out_vertex in incident_edge_list
        for row=1:size(edges,1)
            # if (!(has_edge(tree,(edges[row,:in],edges[row,:out])))) & (indegree(tree, out_vertex) >= 2)
            # TODO: Have AND (if edge.src || edge.dst == out_vertex)
            if (!(has_edge(tree,(edges[row,:in],edges[row,:out]))) .& ((edges[row,:in]==out_vertex) .| (edges[row,:out]==out_vertex)))
                swap_value = edges[row,:weight] - find_weight((vertex,out_vertex), edges)
                # println(find_weight((vertex,out_vertex), edges))
                if swap_value < best_swap_value
                    best_swap_value = swap_value
                    best_swap = [(vertex, out_vertex), (edges[row,:in],edges[row,:out])]
                end
            end
        end
    end
    # println(best_swap)
    # println(best_swap_value)
    rem_edge!(tree, best_swap[1][1], best_swap[1][2])
    add_edge!(tree, best_swap[2])

    return tree

end


function swap_edges_add!(vertex, tree, edges)
    # First find the best two edges to swap
    best_swap_value = Inf
    best_swap = [(1,2),(2,3)]
    function find_weight(edge, edge_df)
        if edge[1] < edge[2]
            return edge_df[(edge_df[:in].==edge[1]) .& (edge_df[:out].==edge[2]),:weight][1]
        elseif edge[1] > edge[2]
            return edge_df[(edge_df[:in].==edge[2]) .& (edge_df[:out].==edge[1]),:weight][1]
        else
            println("ERROR")
        end
    end

    incident_edge_list = map((x,y) -> x==vertex ? y : x, edges[:in],edges[:out])
    incident_edge_list = symdiff(incident_edge_list,outneighbors(tree,vertex))
    # TODO: Fix to ensure edge (1,n+1) is not swapped out
    # println("index list is: ", incident_edge_list)
    for in_vertex in incident_edge_list
        for row=1:size(edges,1)
            # println(edges[row,:])
            # if ((edges[row,:in]!=n+1) & (edges[row,:out]!=n+1))
            # if (has_edge(tree,(edges[row,:in],edges[row,:out])) .& (edges[row,:in]!=n+1) .& (edges[row,:out]!=n+1) .& (edges[row,:in]!=1) .& (edges[row,:out]!=1))
            if (has_edge(tree,(edges[row,:in],edges[row,:out])) .& (edges[row,:in]!=1) .& (edges[row,:out]!=1))
                # println(edges[row,:in]!=n+1)
                # println(edges[row,:out]!=n+1)
                swap_value = edges[row,:weight] - find_weight((vertex,in_vertex), edges)
                if swap_value < best_swap_value
                    # Create a temp graph to test connectedness
                    temp_tree = copy(tree)
                    add_edge!(temp_tree, vertex, in_vertex)
                    set_prop!(temp_tree, vertex, in_vertex, :weight, find_weight((vertex,in_vertex),edges))
                    rem_edge!(temp_tree, edges[row,:in], edges[row,:out])
                    if temp_tree |> is_connected
                        best_swap_value = swap_value
                        best_swap = [(vertex, in_vertex), (edges[row,:in],edges[row,:out])]
                    end
                end
            end
        # end
        end
    end
    # println("Best swap is: ",best_swap)
    add_edge!(tree, best_swap[1])
    set_prop!(tree, best_swap[1][1], best_swap[1][2], :weight, find_weight(best_swap[1], edges))
    rem_edge!(tree, best_swap[2][1], best_swap[2][2])

    return tree

end

function append_new_col!(new_column)
    # columns[τ] = new_column
    # Get objective value of new column
    cost = 0
    # println(distances)
    # println(new_column)
    for v in new_column.vertex
        # println(v, "has cost ", positions[:profit][v])
        cost -= positions[:profit][v]
    end
    for edge in new_column.edge
        # println(edge, "has cost ", distances[edge.src,edge.dst])
        cost += distances[edge.src,edge.dst]
    end
    # println("Cost is: ", cost)


    # Get value for constraint (10a)
    value = OffsetArray{Float64}(3:n)
    for v=3:n
        value[v] = mapreduce(x -> (x.src==v || x.dst == v) ? 1 : 0, +, 0, new_column.edge) - 2* (v in new_column.vertex ? 1 : 0)
    end

    # Get value for constraint (10b)
    weight_value = sum(positions[:weight][v] for v in new_column.vertex if v!=1) - max_weight

    # Add the new column
    # Make the coefficients an array
    coeff_array = []
    for i in [value, weight_value]
        append!(coeff_array,i)
    end
    global convexity_constraint
    if τ==1
        # global convexity_constraint = @constraint(rmp_τ, convexity_constraint, 0 == 1)
        global convexity_constraint = @constraint(rmp_τ, 0 == 1)
        push!(constraint_refs, convexity_constraint)
    end
    @variable(rmp_τ, 0 <= λ_new <= Inf, objective = cost, inconstraints = constraint_refs, coefficients = convert(Array{Float64,1}, [coeff_array;1]))
    setname(λ_new,string("λ[",τ,"]"))
    return cost, coeff_array
end

solve(rmp_τ)
function c1t_problem_solve(w_hat, r_hat)
    c1t_distances = zeros(n+1, n+1)
    c1t_distances[n+1,1:n] = r_hat
    c1t_distances[1:n,n+1] = r_hat
    c1t_distances[1:n,1:n] = w_hat
    distance_graph = SimpleWeightedGraph(n+1)
    edges_df = DataFrame(in = Int64[], out = Int64[], weight = Float64[])
    for i=1:n+1, j=i+1:n+1
        LightGraphs.add_edge!(distance_graph, i, j, c1t_distances[i,j])
        push!(edges_df, [i j c1t_distances[i,j]])
    end
    sort!(edges_df, :weight)
    mst = LightGraphs.kruskal_mst(distance_graph)
    # tree = SimpleWeightedGraph(n+1)
    tree = MetaGraph(n+1)
    for i=1:length(mst)
        add_edge!(tree, mst[i].src, mst[i].dst)
        set_prop!(tree, mst[i].src, mst[i].dst, :weight, distance_graph.weights[mst[i].src,mst[i].dst])
    end
    if !(has_edge(tree,1,n+1))
        add_edge!(tree, 1, n + 1)
        set_prop!(tree, 1, n + 1, :weight, distance_graph.weights[1,n+1])
    else
        for row=1:size(edges_df,1)
            if !(has_edge(tree, (edges_df[row,:in], edges_df[row,:out])))
                add_edge!(tree, edges_df[row,:in], edges_df[row,:out])
                set_prop!(tree, edges_df[row,:in], edges_df[row,:out], :weight, edges_df[row, :weight])
                break
            end
        end
    end
    # println(tree |> edges)
    # println(tree |> edges |> collect)
    iter = 1
    # println("The degree of vertex 1 is: ", indegree(tree,1))
    while indegree(tree,1) != 3
        if indegree(tree, 1) > 3
            swap_edges_remove!(1,tree, edges_df)
        elseif indegree(tree,1) < 3
            swap_edges_add!(1,tree, edges_df)
        end
        if iter > 100
            break
        end
    end
    # println(tree |> edges |> collect)

    feasible_edges = []
    feasible_vertices = Set()
    for edge in tree |> edges
        if ((edge.src!=n+1) .& (edge.dst!=n+1))
            push!(feasible_edges,edge)
            push!(feasible_vertices, edge.src, edge.dst)
        # end
        else
            delete!(feasible_vertices, min(edge.src, edge.dst))
        end
    end
    push!(feasible_vertices,1)
    # return soln(feasible_edges, feasible_vertices)
    return soln(feasible_edges, feasible_vertices), tree
end

function reduced_cost(new_column)
    edge_cost = sum(w_hat[edge.src,edge.dst] for edge in new_column.edge)
    vertex_cost = sum(r_hat[v] for v in new_column.vertex)
    # return edge_cost - vertex_cost - getdual(convexity_constraint)
    return edge_cost - vertex_cost
end

function test_integrality(rmp)
    return mapreduce(x -> abs(getvalue(x))<0.000000001 ? true : false, &, true, Variable.(rmp, 1:(2*n-2)))
end
lower_bound = []
upper_bound = []
objective_array = []
println("Solving the CGP now")
solve(rmp_τ)
global convexity_constraint
for i=1:10
    # println(Variable.(rmp_τ, 1:(2*n)) |> getvalue)
    if ((test_integrality(rmp_τ)) .& (τ!=1))
        v_ub = min(v_ub, getobjectivevalue(rmp_τ))
        println("New upper bound found: ",v_ub)
        append!(upper_bound, v_ub)
        # break
    end
    global w_hat
    global r_hat
    w_hat, r_hat = create_reduced_costs(rmp_τ)
    # println("w_hat and r_hat have been defined")
    new_column,tree = c1t_problem_solve(w_hat, r_hat)
    # println("new_column has been defined")
    columns[i] = (new_column,tree)
    println("Lower bound is: ",v_lb)
    println("objective value is: ", getobjectivevalue(rmp_τ))
    if τ==1
        # v_lb = max(v_lb, getobjectivevalue(rmp_τ))
    else
        println("Reduced cost is: ",getobjectivevalue(rmp_τ) + reduced_cost(new_column) - getdual(convexity_constraint))
        println(reduced_cost(new_column))
        println(getdual(convexity_constraint))
        v_lb = max(v_lb, getobjectivevalue(rmp_τ) + reduced_cost(new_column) - getdual(convexity_constraint))
    end
    # println(v_lb)
    vertex_duals = getdual(vertex_constraints)
    for i=3:n
        α_tilde[i] = vertex_duals[i]
    end
    β_tilde = getdual(knapsack_constraint)
    cost, = append_new_col!(new_column)
    append!(objective_array,cost)
    # println(getdual(convexity_constraint))
    push!(lower_bound, v_lb)
    @objective(rmp_τ, Min, sum( (α_tilde[v] + ϵ_α) * η_ub[v] - (α_tilde[v] - ϵ_α) * η_lb[v] for v=3:n) +
        (β_tilde + ϵ_β)*η_ub_knap - (β_tilde - ϵ_β)*η_lb_knap + sum(objective_array.*Variable.(rmp_τ, (2*n-1):(2*n-1 + τ-1))))
    println("Iteration number: ",τ)
    solve(rmp_τ)
    if τ>2
        if ((columns[τ][1].edge == columns[τ-1][1].edge) & (columns[τ][1].vertex == columns[τ-1][1].vertex))
            # break
        end
    end

    τ += 1
end

# To get variable values
for var in Variable.(rmp_τ, 1:rmp_τ.numCols)
    println(var, ": ", getvalue(var))
end

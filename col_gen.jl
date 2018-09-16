# This file matches all data and creates the relevant dataframe
using CSV
using DataFrames
using CPLEX
using JuMP
# import Graphs
using LightGraphs
using SimpleWeightedGraphs
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

max_weight = Int(round(n/2))

positions = CSV.read("./data/positions_"*string(n)*".csv")

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')./2
end

distances = get_distance(positions)

################################
# Create the RMP_τ^(I^P)
################################
τ = 1
println("Create the RMP")
ϵ_α = 1
ϵ_β = 1
ϵ_μ = 1

α_tilde = fill!(OffsetArray{Float64}(2:n),0)
β_tilde = 0

I_P = Set{Int64}()
v_lb = -Inf
v_ub = Inf
c_1 = Inf

columns = Dict()

rmp_τ = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))

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

solve(rmp_τ)


# α = OffsetArray{Float64}(3:n)
# for i=3:n
#     α[i]=getdual(vertex_constraints[i])
# end
#
# β = getdual(knapsack_constraint)
# Create the CGP

# cgp = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))
#
# @variable(cgp, edge[i=1:n, j=i+1:n], Bin)
# @variable(cgp, vertex[i=1:n], Bin)


# Create the modified distances
# w_hat = distances
# for i=2:n
#     for j=i+1:n
#         w_hat[i,j] -= α_tilde[i] + α_tilde[j]
#     end
# end
#
# r = positions[:profit]
# for v=3:n
#     r[v] = r[v] - 2*α[v] + positions[:weight][v] * β
# end

# @objective(cgp, Min, sum(w_hat[i,j]*edge[i,j] for i=1:n, j=i+1:n) - sum(r[i] * vertex[i] for i=2:n))
#
# # Constraint (6c)
# @constraint(cgp, sum(edge[1,i] for i=2:n) == 2)
# # Constraint (2a)
# @constraint(cgp, sum(edge[i,j] for i=1:n, j=i+1:n) == sum(vertex[i] for i=1:n))
# @constraint(cgp, vertex[1] == 1)
# # Constraints (2b)
# println("Creating the generalised subtour elimination constraints")
# for num in 2:n
#     println("Up to $num cities")
#     for list in combinations(2:n, num)
#         println(list)
#         println(length(list))
#         @constraint(cgp, sum(edge[i,j] for i in list, j in list if i < j) >= 1 + sum(vertex[i] for i in list if i!=1) - length(list))
#     end
# end
#
# solve(cgp)

function create_reduced_costs(rmp)
    α = getdual(vertex_constraints)
    β = getdual(knapsack_constraint)
    w_hat = distances
    for i=3:n
        for j=i+1:n
            w_hat[i,j] -= α[i] + α[j]
        end
    end

    r = positions[:profit]
    r_hat = r
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

## Create the C1T problem
w_hat, r_hat = create_reduced_costs(rmp_τ)
c1t_distances = zeros(n+1, n+1)
c1t_distances[n+1,1:n] = r_hat
c1t_distances[1:n,n+1] = r_hat
c1t_distances[1:n,1:n] = w_hat
distance_graph = SimpleWeightedGraph(n+1)
# distance_graph = LightGraphs.AbstractGraph(n+1)
edges_df = DataFrame(in = Int64[], out = Int64[], weight = Float64[])
for i=1:n+1, j=i+1:n+1
    # Graphs.add_edge!(distance_graph, i, j, c1t_distances[i,j])
    LightGraphs.add_edge!(distance_graph, i, j, c1t_distances[i,j])
    push!(edges_df, [i j c1t_distances[i,j]])
end
sort!(edges_df, :weight)
mst = LightGraphs.kruskal_mst(distance_graph)
# mst = Graphs.kruskal_minimum_spantree(distance_graph)
# tree = Graphs.simple_graph(n+1, is_directed = false)
# mst_distances = copy(c1t_distances)
# for i=1:length(mst)
#     Graphs.add_edge!(tree, mst[i].src, mst[i].dst)
#     if mst[i].src < mst[i].dst
#         mst_distances[mst[i].src, mst[i].dst] = Inf
#     else
#         mst_distances[mst[i].dst, mst[i].src] = Inf
#     end
# end
if !(LightGraphs.Edge(1,n+1) in mst)
    add_edge!(mst, 1, n + 1)
end

function swap_edges_remove(vertex, tree, edges)
    # First find the best two edges to swap
    best_swap_value = Inf
    best_swap = [(1,2),(2,3)]
    function find_weight(edge, edge_df)
        if edge.source < edge.target
            return edge_df[(edge_df[:in].==edge.source) .& (edge_df[:out].==edge.target),:weight]
        elseif edge.source > edge.target
            return edge_df[(edge_df[:in].==edge.target) .& (edge_df[:out].==edge.source),:weight]
        else
            println("ERROR")
        end
    end
    incident_edge_list = Graphs.out_edges(vertex, tree)
    println(incident_edge_list)
    for edge in incident_edge_list
        for row=1:size(edges,1)
            println((edges[row,:in],edges[row,:out]))
            if !((edges[row,:in],edges[row,:out]) in tree.edges)
                swap_value = edges[row,:weight] - find_weight(edge, edges)[1]
                println(swap_value)
                if swap_value < best_swap_value
                    best_swap_value = swap_value
                    best_swap = [(edge.source, edge.target), (edges[row,:in],edges[row,:out])]
                end
            end
        end
    end
    println(best_swap)
    println(best_swap_value)

    temp = sort(edges[((edges[:in].==vertex) .| (edges[:out].==vertex)),:],:weight)

    return tree

end

# while Graphs.in_degree(1, tree) != 3
#     if Graphs.in_degree(1, tree) > 3
#         # Find biggest weight edge leaving 1
#         temp = sort(edges_df[((edges_df[:in].==1) .| (edges_df[:out].==1)),:],:weight)
#     elseif Graphs.in_degree(1, tree) < 3
        temp = sort(edges_df[((edges_df[:in].==1) .| (edges_df[:out].==1)),:],:weight)
#
#     end
# end
#
# min_weight = Inf
# local min_edge
# for i=1:n+1
#     if mst_distances[1,i] <= min_weight
#         min_weight = mst_distances[1,i]
#         min_edge = [1,i]
#     end
# end
# Graphs.add_edge!(tree, min_edge[1], min_edge[2])
#
#
# # new_column = soln(getvalue(edge), getvalue(vertex))
# # Append new column
#
# function append_new_col!(new_column)
#     columns[τ] = new_column
#     # Get objective value of new column
#     cost = 0
#     for i=1:n
#         cost -= positions[:weight][i] * getvalue(vertex[i])
#         for j=i+1:n
#             cost += distances[i,j] * getvalue(edge[i,j])
#         end
#     end
#
#     # Get value for constraint (10a)
#     value = OffsetArray{Float64}(3:n)
#     for k=3:n
#         value[k] = sum(new_column.edge[i,j] for i=1:n, j=1:n if (i<j && (i==k || j==k))) - 2* new_column.vertex[k]
#     end
#
#     # Get value for constraint (10b)
#     weight_value = sum(positions[:weight][vertex] * new_column.vertex[vertex] for vertex=2:n) - max_weight
#
#     # Add the new column
#     # Make the coefficients an array
#     coeff_array = []
#     for i in [value, weight_value]
#         append!(coeff_array,i)
#     end
#     # @variable(rmp_τ, 0 <= λ_new <= Inf, objective = cost, inconstraints = constraint_refs, coefficients = convert(Array{Float64,1},coeff_array))
#     if !(isdefined(:convexity_constraint))
#         @constraint(rmp_τ, convexity_constraint, 0 == 1)
#         push!(constraint_refs, convexity_constraint)
#     end
#     @variable(rmp_τ, 0 <= λ_new <= Inf, objective = cost, inconstraints = constraint_refs, coefficients = convert(Array{Float64,1}, [coeff_array;1]))
#     setname(λ_new,string("λ[",τ,"]"))
# end

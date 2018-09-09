
# This file matches all data and creates the relevant dataframe

using CSV
using DataFrames
using CPLEX
using JuMP
import LightGraphs
using GraphLayout
using Combinatorics

include("helper.jl")

if length(ARGS) >= 1
    println(ARGS[1])
    global n = parse(Int64,ARGS[1])
elseif !isdefined(:n)
    global n = 10
# else
#     global n = 5
end
# max_weight = 7
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

for n in 150:50:1000
    positions = create_problem(n)
    CSV.write("data/positions_"*string(n)*".csv", positions)
end
# 
# distances = get_distance(positions)
#
# tspst = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))
#
# # Add edge variables
# @variable(tspst, edge[i=1:n,j=i+1:n], Bin)
# @variable(tspst, vertex[i=1:n], Bin)
#
# @objective(tspst, Min, sum(distances[i,j].*edge[i,j] for i=1:n, j=i+1:n) - sum(positions[i,4].*vertex[i] for i=1:n))
#
# # Include the first node
# @constraint(tspst, vertex[1] == 1)
#
# # Only choose max_weight cities
# @constraint(tspst, sum(positions[i,5]*vertex[i] for i=1:n) <= max_weight)
#
# # GSEC
# println("Creating the generalised subtour elimination constraints")
# for num in 3:n
#     println("Up to $num cities")
#     for list in combinations(1:n, num)
#         for k in list
#             if !(1 in list)
#                 @constraint(tspst, sum(edge[i,j] for i in list, j in list if i < j) - sum(vertex[i] for i in list if i!=1) <= - vertex[k])
#             end
#         end
#     end
# end
#
#
# # Vertex degree restrictions
# for i=1:n
#     @constraint(tspst, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2*vertex[i])
# end
#
# println("Solving the model now")
# solve(tspst)
#
# x = zeros(n,n)
# edge_values = getvalue(edge)
# for i=1:n, j=i+1:n
#     x[i,j] = edge_values[i,j]
# end
# x = x + x'
#
# draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*".svg")
#
#
# ################################
# # Create the RMP_τ
# ################################
#
# rmp_τ = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))

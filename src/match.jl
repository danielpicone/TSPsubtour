# This file matches all data and creates the relevant dataframe
# module match

using CSV
# using DataFrames
using CPLEX
using JuMP
import LightGraphs
# using Combinatorics

positions = CSV.read("../data/a03. Position.csv", header = 3)[:,2:end]

times = CSV.read("../data/a02. Time.csv", header = 3)[:,2:end]

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')./2
end

distances = get_distance(positions)

n = size(positions)[1]
tsptw = JuMP.Model(solver = CplexSolver())

# Add edge variables
@variable(tsptw, edge[i=1:n,j=i+1:n], Bin)
@variable(tsptw, vertex[i=1:n], Bin)
# @variable(tsptw, 0 <= time <= maximum(convert(Array{Float64,2},times[:,2:3])))
@variable(tsptw, 0 <= time <= Inf)

@objective(tsptw, Min, sum(distances[i,j].*edge[i,j] for i=1:n, j=i+1:n) - sum(positions[i,4].*vertex[i] for i=1:n))

# Include the first node
@constraint(tsptw, vertex[1] == 1)

# Leave every node once
# for i=1:n
#     @constraint(tsptw, sum(edge[i,j] for j=1:n if j!=i) == 1)
# end

# GSEC
for num in 3:n
    # for list in Combinatorics.permutations(1:n, num)
    for list in combinations(1:n, num)
        for k in list
            if !(1 in list)
                println(list)
                # @constraint(tsptw, sum(edge[i,j] for i in list, j in list if i != j) - 1*sum(vertex[i] for i in list) <= - 1*vertex[k])
                @constraint(tsptw, sum(edge[i,j] for i in list, j in list if i < j) - 1*sum(vertex[i] for i in list if i!=1) <=0)
            end
        end
    end
end

# Only choose n cities
@constraint(tsptw, sum(vertex[i] for i=1:n) <= 3)

# Vertex degree restrictions
for i=1:n
    @constraint(tsptw, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2*vertex[i])
end

solve(tsptw)

x = zeros(n,n)
edge_values = getvalue(edge)
for i=1:n, j=i+1:n
    x[i,j] = edge_values[i,j]
end
x = x + x'

draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="graph.svg")

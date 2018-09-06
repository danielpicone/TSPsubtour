
# This file matches all data and creates the relevant dataframe

using CSV
# using DataFrames
using CPLEX
using JuMP
import LightGraphs
using GraphLayout
using Combinatorics

positions = CSV.read("./data/a03. Position.csv", header = 3)[:,2:end]

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')./2
end

distances = get_distance(positions)

global n = size(positions)[1]
tsptw = JuMP.Model(solver = CplexSolver())

# Add edge variables
@variable(tsptw, edge[i=1:n,j=i+1:n], Bin)
@variable(tsptw, vertex[i=1:n], Bin)

@objective(tsptw, Min, sum(distances[i,j].*edge[i,j] for i=1:n, j=i+1:n) - sum(positions[i,4].*vertex[i] for i=1:n))

# Include the first node
@constraint(tsptw, vertex[1] == 1)

# GSEC
for num in 3:n
    # for list in Combinatorics.permutations(1:n, num)
    for list in combinations(1:n, num)
        for k in list
            if !(1 in list)
                # @constraint(tsptw, sum(edge[i,j] for i in list, j in list if i != j) - 1*sum(vertex[i] for i in list) <= - 1*vertex[k])
                @constraint(tsptw, sum(edge[i,j] for i in list, j in list if i < j) - 1*sum(vertex[i] for i in list if i!=1) <=0)
            end
        end
    end
end

# Only choose n cities
@constraint(tsptw, sum(vertex[i] for i=1:n) <= 4)

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

# draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="graph.svg")
# draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph.svg")

# Create an adjacency matrix from a tour
function get_adj_mat(tour)
    if tour[1]!=1
        println("Error: The tour must start at 1")
    end
    # adj = zeros(n,n)
    adj = spzeros(n,n)
    for i=1:(length(tour)-1)
        adj[tour[i],tour[i+1]] = 1
    end
    adj[1,tour[end]] = 1
    return adj
end
# Create a function which returns the objective cost of a tour
function get_cost(tour, distances::Array{Float64,2})
    cost = 0
    for city in tour
        cost += positions[:profit][city]
    end
    return -cost + get_adj_mat(tour).*distances
end

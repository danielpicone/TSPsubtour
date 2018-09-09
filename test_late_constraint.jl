
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
    const global n = parse(Int64,ARGS[1])
elseif !isdefined(:n)
    const global n = 10
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

# for l=1:10
    # println("Test $l")
function solve_late_con(n)
    positions = create_csv(n)

    distances = get_distance(positions)

    global tspst = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 0))

    # Add edge variables
    @variable(tspst, edge[i=1:n,j=i+1:n], Bin)
    @variable(tspst, vertex[i=1:n], Bin)

    @objective(tspst, Min, sum(distances[i,j].*edge[i,j] for i=1:n, j=i+1:n) - sum(positions[i,4].*vertex[i] for i=1:n))

    # Include the first node
    @constraint(tspst, vertex[1] == 1)

    # Only choose max_num_cities cities
    @constraint(tspst, sum(positions[i,5]*vertex[i] for i=1:n) <= max_weight)

    # Vertex degree restrictions
    for i=1:n
        @constraint(tspst, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2*vertex[i])
    end

    println("Solving the model now")
    solve(tspst)

    # draw_layout_adj(get_adj_mat(edge), convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")
    # draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")

    # Determine if multiple subtours exist

    function get_subtour(edge_solution)
        list_of_edges = []
        for i=1:n, j=i+1:n
            if abs(edge_solution[i,j]-1) < 0.000001
                push!(list_of_edges,(i,j))
            end
        end

        function get_next_city(city, edge)
            if edge[1]==city
                return edge[2]
            else
                return edge[1]
            end
        end

        tour = Array{Int64,1}(0)
        next_city = 1
        push!(tour, next_city)
        while length(list_of_edges) > 0
            index = 1
            for pair in list_of_edges
                if next_city in list_of_edges[index]
                    next_city = get_next_city(next_city, list_of_edges[index])
                    push!(tour, next_city)
                    deleteat!(list_of_edges, index)
                    break
                end
                index += 1
            end
        end
        return tour[1:end-1]
    end

    while !is_valid_solution(edge)
        subtours = get_subtours(edge)
        println("The number of subtours is: ", length(subtours))
        if length(subtours)!=1
            for tour in subtours
                if !(1 in tour)
                    add_gsec!(tour)
                end
            end
        end
        solve(tspst)
        # println(MathProgBase.numlinconstr(tspst))
    end

    # draw_layout_adj(get_adj_mat(edge), convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late_solved.svg")


end

for l=1:10
    println(l)
    # positions = CSV.read("./data/positions_"*string(n)*".csv")
    @time solve_late_con(n)
end

# This file creates the benchmarks scenarios

using CSV
using DataFrames
using CPLEX
using JuMP
import LightGraphs
using GraphLayout
using Combinatorics

include("helper.jl")

# if length(ARGS) >= 1
#     println(ARGS[1])
#     global n = parse(Int64,ARGS[1])
# elseif !isdefined(:n)
#     global n = 10
# end

# positions = CSV.read("./data/positions_"*string(n)*"_v1.csv")

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')
end

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

# function add_gsec_tsp!(subtour)
#     # println(subtour)
#     for k in subtour
#         @constraint(tsp, sum(edge[i,j] for i in subtour,j in subtour if i < j) <= length(subtour)-1)
#     end
# end
# distances = get_distance(positions)
#
# tsp = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))
#
# # Add edge variables
# @variable(tsp, edge[i=1:n,j=i+1:n], Bin)
#
# @objective(tsp, Min, sum(distances[i,j]*edge[i,j] for i=1:n, j=i+1:n))
#
# # Vertex degree restrictions
# for i=1:n
#     @constraint(tsp, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2)
# end
#
# println("Solving the model now")
# solve(tsp)
#
# # draw_layout_adj(get_adj_mat(edge), convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")
# # draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")
#
# # Determine if multiple subtours exist
#
# start_time = time()
# while !is_valid_solution(edge)
#     subtours = get_subtours(edge)
#     if length(subtours)!=1
#         for tour in subtours
#             if !(1 in tour)
#                 add_gsec_tsp!(tour)
#             end
#         end
#     end
#     solve(tsp)
#     # println(MathProgBase.numlinconstr(tsp))
# end
# end_time = time()
#
# println("Time taken was: ", end_time - start_time, " seconds")
# draw_layout_adj(get_adj_mat(edge), convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late_solved.svg")

function solve_tsp(n, v)
    global len = n
    positions = CSV.read("./data/positions_"*string(n)*"_v"*string(v)*".csv")
    distances = get_distance(positions)

    function add_gsec_tsp!(subtour)
        # println(subtour)
        for k in subtour
            @constraint(tsp, sum(edge[i,j] for i in subtour,j in subtour if i < j) <= length(subtour)-1)
        end
    end
    tsp = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 0, CPXPARAM_MIP_Display = 2))
    # Add edge variables
    @variable(tsp, edge[i=1:n,j=i+1:n], Bin)
    @objective(tsp, Min, sum(distances[i,j]*edge[i,j] for i=1:n, j=i+1:n))
    # Vertex degree restrictions
    for i=1:n
        @constraint(tsp, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2)
    end
    solve(tsp)
    while !is_valid_solution(edge)
        subtours = get_subtours(edge)
        if length(subtours)!=1
            for tour in subtours
                if !(1 in tour)
                    add_gsec_tsp!(tour)
                end
            end
        end
        solve(tsp)
    end
    solution_path = []
    valx = getvalue(edge)
    for i=1:n,j=i+1:n
        if valx[i,j]==1
            push!(solution_path, (i,j))
        end
    end
    return getobjectivevalue(tsp), solution_path
end

df = DataFrame(num = Int64[], v = Int64[], value = Float64[], solution_path = Array[], time = Float64[])
for k=10:5:200
    println("Up to number: ",k)
    avg_time = []
    for v=1:10
        start_time = time()
        value, path = solve_tsp(k, v)
        end_time = time()
        time_taken = end_time - start_time
        # df[:n] = k
        # df[:version] = v
        # df[:solution_value] = value
        # df[:solution_path] = [path]
        push!(df, [k, v, value, path, time_taken])
        push!(avg_time, time_taken)
    end
    println(mean(avg_time))
end

CSV.write("./data/solutions.csv", df)

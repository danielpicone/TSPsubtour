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
end
max_num_cities = 7

positions = CSV.read("./data/positions_"*string(n)*".csv")

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')./2
end

function create_problem(num_cities)
    df = DataFrame([1:num_cities rand(-num_cities:num_cities,num_cities,2) rand(1:num_cities,num_cities)])
    df[1,2:end] = 0
    rename!(df, :x1 => :city, :x2 => :xcoord, :x3 => :ycoord, :x4 => :profit)
    return df
end


distances = get_distance(positions)

tspst = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 1, CPXPARAM_MIP_Display = 2))

# Add edge variables
@variable(tspst, edge[i=1:n,j=i+1:n], Bin)
@variable(tspst, vertex[i=1:n], Bin)

@objective(tspst, Min, sum(distances[i,j].*edge[i,j] for i=1:n, j=i+1:n) - sum(positions[i,4].*vertex[i] for i=1:n))

# Include the first node
@constraint(tspst, vertex[1] == 1)

# Only choose max_num_cities cities
@constraint(tspst, sum(vertex[i] for i=1:n) <= max_num_cities)



# Vertex degree restrictions
for i=1:n
    @constraint(tspst, sum(edge[i,j] for j=i+1:n if i!=j) + sum(edge[j,i] for j=1:n if j<i) == 2*vertex[i])
end

println("Solving the model now")
solve(tspst)

draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")

x = zeros(n,n)
edge_values = getvalue(edge)
for i=1:n, j=i+1:n
    x[i,j] = edge_values[i,j]
end
x = x + x'

# draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late.svg")

# Determine if multiple subtours exist

function get_subtour(edge_solution)
    list_of_edges = []
    for i=1:n, j=i+1:n
        if(edge_solution[i,j])==1
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

function partition_edges(edge_solution::JuMP.JuMPDict{JuMP.Variable,2})
    return partition_edges(getvalue(edge_solution))
end

function partition_edges(edge_solution)
    list_of_edges = []
    for i=1:n, j=i+1:n
        if(edge_solution[i,j])==1
            push!(list_of_edges,(i,j))
        end
    end
    tours = Array{Array{Tuple{Int64, Int64},1}}(n)

    function test_if_edge_in_set(vertex_set, edge)
        if edge[1] in vertex_set || edge[2] in vertex_set
            return true
        else
            return false
        end
    end

    tours = Dict()
    num_tours = 1
    tours[1] = Dict("edges" => [list_of_edges[1]], "set" => Set(list_of_edges[1]))
    deleteat!(list_of_edges,1)
    for edge in list_of_edges
        is_in_existing_tour = false
        for index=1:num_tours
            is_in_existing_tour = test_if_edge_in_set(tours[index]["set"], edge)
            if is_in_existing_tour
                push!(tours[index]["edges"],edge)
                push!(tours[index]["set"], edge[1], edge[2])
                break
            end
        end
        if !is_in_existing_tour
            num_tours += 1
            tours[num_tours] = Dict("edges" => [edge], "set" => Set(edge))
        end
    end
    return tours

end

println("hey")
subtours_dict = edge |> partition_edges |> consolidate_sets
subtours = Array{Array{Int64,1},1}(length(subtours_dict))
for (index, key) in enumerate(keys(subtours_dict))
    subtours[index] = create_tour(subtours_dict[key]["edges"])
end

function add_gsec!(subtour)
    for k in subtour
        @constraint(tspst, sum(edge[i,j] for i in subtour,j in subtour if i < j) <= sum(vertex[i] for i in subtour if i != k))
    end
end

function is_valid_solution(edge_solution)
    return edge_solution |> getvalue |> is_valid_solution
end

function is_valid_solution(edge_solution)
    subtours_dict = edge |> partition_edges |> consolidate_sets
    subtours = Array{Array{Int64,1},1}(length(subtours_dict))
    for (index, key) in enumerate(keys(subtours_dict))
        subtours[index] = create_tour(subtours_dict[key]["edges"])
    end
    if subtours |> length > 1
        return false
    elseif subtours |> length == 0
        return true
    else
        println("ERROR: The length of subtours is ", length(subtour))
    end
end

for i=1:100
    solve(tspst)
    for tour in subtours
        if !(1 in tour)
            add_gsec!(tour)
        end
    end
end

draw_layout_adj(x, convert(Array{Float64},positions[:xcoord]), convert(Array{Float64},positions[:ycoord]), filename="./graphs/graph_"*string(n)*"_late_solved.svg")

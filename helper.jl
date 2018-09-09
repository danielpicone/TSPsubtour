# This is a module which provides helpful functions for getting and creating solutions
# module helper
# using JuMP
# export get_subtour, get_adj_mat, get_cost

function create_csv(num_cities)
    df = DataFrame([1:num_cities rand(-num_cities:num_cities,num_cities,2) rand(1:num_cities,num_cities) rand(1:Int(round(num_cities/2)), num_cities)])
    df[1,2:end] = 0
    rename!(df, :x1 => :city, :x2 => :xcoord, :x3 => :ycoord, :x4 => :profit, :x5 => :weight)
    return df
end

function get_subtour(edge_solution::JuMP.JuMPDict{JuMP.Variable,2})
    return get_subtour(getvalue(edge_solution))
end

function get_subtour(edge_solution)
    list_of_edges = []
    for i=1:n, j=i+1:n
        if(edge_solution[i,j])==1
            push!(list_of_edges,(i,j))
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

function get_next_city(city, edge)
    if edge[1]==city
        return edge[2]
    else
        return edge[1]
    end
end

# TODO: Fix this *think* the problem is that if the next city isn't in the first
# index it does not terminate
function create_tour(list_of_edges)
    tour = Array{Int64,1}(0)
    # next_city = 1
    next_city = list_of_edges[1][1]
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
# Create an adjacency matrix from a tour
function get_adj_mat(tour::JuMP.JuMPDict{JuMP.Variable,2})
    x = zeros(n,n)
    edge_values = getvalue(tour)
    for i=1:n, j=i+1:n
        x[i,j] = edge_values[i,j]
    end
    x = x + x'
    return x
end

function get_adj_mat(tour)
    if tour[1]!=1
        println("Error: The tour must start at 1")
    end
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
    return -cost + sum(get_adj_mat(tour).*distances)
end

function consolidate_sets(dict)
    len = length(dict)
    for i=1:len
        for j=i+1:len
            if haskey(dict,i) & haskey(dict,j)
                if length(intersect(dict[i]["set"], dict[j]["set"])) >= 1
                    append!(dict[i]["edges"], dict[j]["edges"])
                    dict[i]["set"] = union(dict[i]["set"],dict[j]["set"])
                    delete!(dict, j)
                end
            end
        end
    end
    return dict
end

function get_subtours(edge)
    subtours_dict = edge |> partition_edges |> consolidate_sets
    subtours = Array{Array{Int64,1},1}(length(subtours_dict))
    for (index, key) in enumerate(keys(subtours_dict))
        subtours[index] = create_tour(subtours_dict[key]["edges"])
    end
    return subtours
end

function add_gsec!(subtour)
    for k in subtour
        @constraint(tspst, sum(edge[i,j] for i in subtour,j in subtour if i < j) <= sum(vertex[i] for i in subtour if i != k))
    end
end

function partition_edges(edge_solution::JuMP.JuMPDict{JuMP.Variable,2})
    return partition_edges(getvalue(edge_solution))
end

function partition_edges(edge_solution)
    list_of_edges = []
    for i=1:n, j=i+1:n
        # if(edge_solution[i,j])==1
        if abs(edge_solution[i,j]-1) < 0.000001
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

function partition_edges(edge_solution::JuMP.JuMPDict{JuMP.Variable,2})
    return partition_edges(getvalue(edge_solution))
end

function partition_edges(edge_solution)
    list_of_edges = []
    for i=1:n, j=i+1:n
        # if(edge_solution[i,j])==1
        if abs(edge_solution[i,j]-1) < 0.000001
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

function consolidate_sets(dict)
    len = length(dict)
    consolidated_flag = true
    while consolidated_flag
        consolidated_flag = false
        for i=1:len
            for j=i+1:len
                if haskey(dict,i) & haskey(dict,j)
                    if length(intersect(dict[i]["set"], dict[j]["set"])) >= 1
                        append!(dict[i]["edges"], dict[j]["edges"])
                        dict[i]["set"] = union(dict[i]["set"],dict[j]["set"])
                        delete!(dict, j)
                        consolidated_flag = true
                    end
                end
            end
        end
    end
    return dict
end

function is_valid_solution(edge_solution::JuMP.JuMPDict{JuMP.Variable,2})
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
    elseif subtours |> length == 1
        return true
    else
        println("ERROR: The length of subtours is ", length(subtours))
    end
end

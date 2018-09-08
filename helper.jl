# This is a module which provides helpful functions for getting and creating solutions
# module helper
# using JuMP
# export get_subtour, get_adj_mat, get_cost

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
    for i=1:length(dict)
        for j=i+1:length(dict)
            if length(intersect(dict[i]["set"], dict[j]["set"])) >= 1
                append!(dict[i]["edges"], dict[j]["edges"])
                dict[i]["set"] = union(dict[i]["set"],dict[j]["set"])
                delete!(dict, j)
            end
        end
    end
    return dict
end
# end

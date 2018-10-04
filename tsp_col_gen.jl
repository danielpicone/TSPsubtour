
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
using GLPKMathProgInterface

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

positions = CSV.read("./data/positions_"*string(n)*"_v2.csv")
positions[:weight] = ones(n)
positions[:profit] = positions[:profit]*3

function get_distance(coord)
    n = size(coord)[1]
    D = zeros(n,n)
    for i=1:n, j=i:n
        D[i,j] = sqrt((coord[i,2]-coord[j,2])^2 + (coord[i,3] - coord[j,3])^2)
    end
    return (D+D')
end

distances = get_distance(positions)

if n==5
    distances = [0 7 2 1 5;
                 7 0 3 6 8;
                 2 3 0 4 2;
                 1 6 4 0 9;
                 5 8 2 9 0]
    cons_mat = [2 2 2 2 2 2 2;
                2 2 2 1 1 2 3;
                2 3 2 3 2 3 1;
                2 2 3 3 3 1 1;
                2 1 1 1 2 2 3]
    obj_array = [28 25 21 19 22 18 28]
end
################################
# Create the RMP_τ^(I^P)
################################
global τ = 1
println("Create the RMP")

# Get cost of inital tour
initial_cost = sum(distances[i,i+1] for i=1:n-1) + distances[1,n]
columns = Dict()

rmp_τ = JuMP.Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 0, CPXPARAM_Preprocessing_Dual=-1, CPXPARAM_MIP_Display = 3))
# rmp_τ = Model(solver = GLPKSolverLP(method=:Exact,msg_lev=3))

@variable(rmp_τ, λ >= 0)
# @variable(rmp_τ, λ[1:7] >= 0)

# Add the objective function
@objective(rmp_τ, Min, initial_cost * λ)
# @objective(rmp_τ, Min, sum(obj_array[i] * λ[i] for i=1:7))

@constraint(rmp_τ, vertex_constraints[i=1:n], 2*λ == 2)
# @constraint(rmp_τ, vertex_constraints[i=1:n], sum(cons_mat[i,j]*λ[j] for j=1:7) == 2)

function create_reduced_costs(rmp)
    α = getdual(vertex_constraints)
    d_hat = copy(convert(Array{Float64,2},distances))
    for i=1:n
        for j=i+1:n
            d_hat[i,j] += -α[i] - α[j]
            d_hat[j,i] = d_hat[i,j]
        end
    end

    return d_hat
end

function find_one_tree(d_hat, distances; num_trees = 1)
    graph = MetaGraph(n-1)
    for i=1:n-1, j=i+1:n-1
        add_edge!(graph, i, j)
        set_prop!(graph, i, j, :weight, d_hat[i+1,j+1])
    end
    tree = kruskal_mst(graph)
    num_v = zeros(n)
    num_v[1] = 2
    # Find number of vertices
    for edge in tree
        num_v[edge.src+1] += 1
        num_v[edge.dst+1] += 1
    end

    cost = 0
    reduced_cost = 0
    for edge in tree
        if edge.src < edge.dst
            cost += distances[edge.src+1,edge.dst+1]
            reduced_cost += d_hat[edge.src+1,edge.dst+1]
        elseif edge.src > edge.dst
            cost += distances[edge.dst+1,edge.src+1]
            reduced_cost += d_hat[edge.dst+1,edge.src+1]
        else
            println("ERROR")
        end
    end
    new_edges = [d_hat[1,2:end] collect(2:n)]'
    new_edges = sortcols(new_edges)
    new_dual = getdual(vertex_constraints)
    # ind1, ind2 = new_edges[2,1:2]
    inds = Int64.(new_edges[2,:])
    # cost += distances[1,Int64(ind1)] + distances[1,Int64(ind2)]
    tree_cost = []
    num_v_array = []
    reduced_cost = []
    for i=1:num_trees
        temp_cost = cost + distances[1,inds[1]] + distances[1,inds[i+1]]
        temp_num_v = copy(num_v)
        temp_num_v[inds[1]] += 1
        temp_num_v[inds[i+1]] += 1
        temp_reduced_cost = temp_cost- sum(temp_num_v[j]*new_dual[j] for j=1:n)
        if temp_reduced_cost > 0
            println("this occurs")
            break
        end
        push!(tree_cost, temp_cost)
        push!(num_v_array, temp_num_v)
        # push!(num_v_array, copy(num_v))
        # println(tree_cost[i])
        # num_v_array[i][inds[1]] += 1
        # num_v_array[i][inds[i+1]] += 1
        push!(reduced_cost, temp_reduced_cost)
    end
    # println("Current cost is: ", cost)
    # println("Cost vector is: ", new_edges)


    return tree_cost, reduced_cost, tree, num_v_array
    # return cost, reduced_cost, tree, num_v
end

function append_new_col!(rmp_τ)
    d_hat = create_reduced_costs(rmp_τ)
    # show(STDOUT,"text/plain", d_hat)
    # println()
    # new_cost, reduced_cost, tree, num_v = find_one_tree(d_hat, distances)
    num_new_cols = 1
    new_costs, reduced_cost, tree, num_v_array = find_one_tree(d_hat, distances; num_trees = num_new_cols)
    # new_cost = find_cost(tree, d_hat, distances)
    # @variable(rmp_τ, 0 <= λ_new <= Inf, objective = new_cost, inconstraints = vertex_constraints, coefficients = num_v)
    # for i=1:num_new_cols
    for (index, new_col_cost) in enumerate(new_costs)
        @variable(rmp_τ, 0 <= λ_new <= Inf, objective = new_col_cost, inconstraints = vertex_constraints, coefficients = num_v_array[index])
        setname(λ_new,string("λ[",τ,"]"))
    end
    # return new_cost, reduced_cost, num_v, tree
    return new_costs, reduced_cost, num_v_array, tree
end

function test_integrality(rmp)
    return mapreduce(x -> abs(getvalue(x))<0.000000001 ? true : false, &, true, Variable.(rmp, 1:(2*n-2)))
end
function print_solution(rmp_τ)
    if round.(getvalue(λ),5)!=1
        for var in Variable.(rmp_τ, 1:rmp_τ.numCols)
            println(var, ": ", getvalue(var))
        end
    end
end
function in_box(dual, center, box)
    return dual < center-box ? center-box : (dual > center+box ? center+box : dual)
end
lower_bound = []
upper_bound = []
objective_array = []
duals = []
println("Solving the CGP now")

# d_hat = create_reduced_costs(rmp_τ)
# tree, num_v = find_one_tree(d_hat)
# new_cost = find_cost(tree, distances)
flag = true
start_time = time()
v_ub = Inf
v_lb = -Inf
# gap = (v_ub-v_lb)/v_lb
gap = Inf
ϵ = 0.0001
# for i=1:4
temp = 0
# box = 1000000000000
center = zeros(n)
# new_dual = optdual
while gap > ϵ
    solve(rmp_τ)
    append!(upper_bound, getobjectivevalue(rmp_τ))
    v_ub = getobjectivevalue(rmp_τ)
    # if τ>2
    new_dual = getdual(vertex_constraints)
        # new_dual = in_box.(new_dual, center, box)
        # center = new_dual
    # end
    push!(duals, new_dual)
    # if τ>2
    #     println("Dual changed by: ", norm(duals[τ]-duals[τ-1]))
    # end
    # println(norm(new_dual,1))
    # print_solution(rmp_τ)
    new_costs, reduced_cost, num_v_array, tree = append_new_col!(rmp_τ)
    temp = copy(reduced_cost)
    # println(new_cost - sum(num_v[i]*new_dual[i] for i=1:n))
    # reduced_cost = new_costs[1] - sum(num_v_array[1][j]*new_dual[j] for j=1:n)
    # println("Last objective value was: ", getobjectivevalue(rmp_τ))
    append!(lower_bound, v_ub+reduced_cost[1])
    v_lb = max(v_lb,v_ub + reduced_cost[1])
    println(gap)
    gap = (v_ub-v_lb)/v_lb
    columns[τ] = tree
    τ+=1
end

solve(rmp_τ)
new_dual = getdual(vertex_constraints)
push!(duals, new_dual)
# append!(lower_bound, v_ub+reduced_cost[1])
# append!(upper_bound, getobjectivevalue(rmp_τ))
end_time = time()
# for var in Variable.(rmp_τ, 1:rmp_τ.numCols)
#     println(var, ": ", getvalue(var))
# end
println("Completed ",τ," iterations in: ", end_time - start_time, " seconds")
relaxed_solution_vars = Variable.(rmp_τ,1:rmp_τ.numCols)[rmp_τ.colVal .!= 0]
for sol in relaxed_solution_vars
    println(sol, " has a value of ", getvalue(sol))
end

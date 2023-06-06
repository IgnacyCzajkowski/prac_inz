using Graphs
using GraphPlot
using Setfield
using Plots

global const inf_prob = 0.5
global const death_prob = 0.0

struct Network
    graph::Graph
    network_state::Vector{Int}
    source_idx::Int
end    

function initialize_info_source(v::Vector)
    n = length(v)
    idx = rand(1:n)
    v[idx] = 1
    return idx #tu bylo dodane
end    

function generate_network(nodes::Int, prob::Float16)
    G::Graph = erdos_renyi(nodes, prob)
    network_state = zeros(nv(G))
    source_idx::Int = initialize_info_source(network_state) #Tu dodane
    N = Network(G, network_state, source_idx)
    return N
end 

function generate_network(nodes::Int, base_nodes::Int, edges::Int)
    G = barabasi_albert(nodes, base_nodes, edges)
    network_state = zeros(nv(G))
    source_idx = initialize_info_source(network_state)#
    N = Network(G, network_state, source_idx)
    return N
end    

function show_network(N::Network, population::Int)
    #println(all_neighbors(N.graph, 1))
    gplot(N.graph, nodelabel=1:population)
    #println("Rysuje")
    #println(population)
end

function interact_witch_closest(N::Network, indx::Int)
    neighbors = all_neighbors(N.graph, indx)
    if N.network_state[indx] == 0
        for i in neighbors
                if N.network_state[i] == 1
                    if rand() < inf_prob
                        return 1
                    end
            end    
        end
    elseif (N.network_state[indx] == 1)
        if rand() < death_prob
            return 2
        end
    end
    return N.network_state[indx]
end    

function get_next_step(N::Network)
    v_next::Vector = zeros(nv(N.graph))
    for indx in 1:length(N.network_state)
        v_next[indx] = N.network_state[indx]
    end    
    for indx in 1:length(N.network_state)
        v_next[indx] = interact_witch_closest(N, indx)
    end    
    return v_next
end 

function get_dist_data(N::Network, ds)
    k::Int = 0
    for i in 1:length(N.network_state)
        if N.network_state[i] == 1
            if ds.dists[i] > k
                k = ds.dists[i]
            end 
        end 
    end 
    return k             
end

function algorithm(N::Network, time_steps::Int)
    dist = Vector{Int}()
    ds = desopo_pape_shortest_paths(N.graph, N.source_idx)
    
    for t in 1:time_steps
        d = get_dist_data(N, ds)
        push!(dist, d)
        N = @set N.network_state = get_next_step(N)
        println("Pokolenie: ", t)
        
        #println(d, "  94")
        #push!(dist, d)
        #println(N.network_state)
    end
    println("Maksymalna odlegosc od zrodla: ", maximum(ds.dists))
    time = collect(1:time_steps)
    #println(time)
    #println(dist)
    scatter(time, dist)
end        


function main()

    params = [1000, 100, 45]
    population::Int = params[1]

    if length(params) == 2
        nodes::Int = params[1]
        prob::Float16  = params[2]
        N::Network = generate_network(nodes, prob) #
    elseif length(params) == 3
        nodes = params[1]
        base_nodes::Int = params[2]
        edges::Int = params[3]
        N = generate_network(nodes, base_nodes, edges) #
        #show_network(N, population)
    end 
    #println(all_neighbors(N.graph, N.source_idx), " Podobno sasiedzi")
    #println(N.network_state) 
    gplot(N.graph, nodelabel=1:population) 
    algorithm(N, 25)     
end 

main()


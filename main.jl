using Graphs
using GraphPlot
using Setfield
using Plots
using Statistics
using StatsBase
using LinearAlgebra

#Struktura implementująca sieć
struct Network
    graph::Graph
    network_state::Vector{Int}
    source_idx::Int
end

#Struktura implementująca obserwatora
mutable struct Observer
    idx::Int
    t::Int
end

#Funkcja inicjalizująca początkowy stan sieci
function initialize_info_source(v::Vector)
    n = length(v)
    idx = rand(1:n)
    v[idx] = 1
    return idx 
end

#Funkcja tworząca sieć typu E-R o podanych parametrach i inicjalizująca stan początkowy sieci
function generate_network(nodes::Int, prob::Float16)
    G::Graph = erdos_renyi(nodes, prob)
    network_state = zeros(nv(G))
    for idx in 1:length(network_state)
        #Zapewnienie braku odizolowanych węzłów
        if length(all_neighbors(G, idx)) == 0
            rand_idx = idx
            while   rand_idx == idx
                rand_idx = rand(1:length(network_state))
            end     
            add_edge!(G, idx, rand_idx)
        end    
    end    
    source_idx::Int = initialize_info_source(network_state) 
    N = Network(G, network_state, source_idx)
    return N
end

#Funkcja tworząca sieć typu B-A o podanych parametrach i inicjalizująca stan początkowy sieci
function generate_network(nodes::Int, base_nodes::Int, edges::Int)
    G = barabasi_albert(nodes, base_nodes, edges, complete = true)
    network_state = zeros(nv(G))
    source_idx = initialize_info_source(network_state)
    N = Network(G, network_state, source_idx)
    return N
end

#Funkcja do aktualizacji stanu węzła o indeksie indx w modelu SI
function  interact_witch_closest(N::Network, indx::Int, inf_prob_loc::Float64)
    neighbors = all_neighbors(N.graph, indx)
    k_interact = 0
    if N.network_state[indx] == 0
        for i in neighbors
            if N.network_state[i] == 1
                k_interact = k_interact + 1
            end
        end
        if rand() < 1 - (1 - inf_prob_loc) ^ k_interact
            return 1
        end
    end
    return N.network_state[indx]
end

#Funkcja do aktualizacji stanu węzła o indeksie indx w modelu FSIR
function  interact_witch_closest_fsir(N::Network, indx::Int, inf_prob_loc::Float64, gamma::Float64)
    neighbors = all_neighbors(N.graph, indx)
    k_interact = 0
    if N.network_state[indx] == 0
        for i in neighbors
            if N.network_state[i] == 1
                k_interact = k_interact + 1
            end
        end
        if rand() < 1 - (1 - inf_prob_loc / k_interact^gamma) ^ k_interact
            return 1
        end
    end
    return N.network_state[indx]
end

#Funkcja wyznaczająca stany węzłów w następnym kroku symulacji
function get_next_step(N::Network, inf_prob_loc::Float64, gamma::Float64)
    v_next::Vector = zeros(nv(N.graph))
    for indx in 1:length(N.network_state) 
        v_next[indx] = N.network_state[indx]
    end
    for indx in 1:length(N.network_state)
        v_next[indx] = interact_witch_closest_fsir(N, indx, inf_prob_loc, gamma)
    end    
    return v_next
end

#Funkcja losująca i inicjalizująca wektor obserwatorów
function getObservers(N::Network, l::Int)
    obs = Vector{Observer}()
    a = sample(1:length(N.network_state), l, replace=false)
    for idx in a
        if idx == N.source_idx
            o = Observer(idx, 0)  
        else    
            o = Observer(idx, Int(floatmax(Float16)))  
        end
        push!(obs, o)
    end
    return obs
end

#Funkcja służąca aktualizowanie stanu obserwatorów w trakcie symulacji
function actuateObservers(N_new::Network, N_old_state::Vector{Int}, obs::Vector{Observer}, time::Int) 
    for point in obs
        if point.t == Int(floatmax(Float16))  
            if N_new.network_state[point.idx] == 1 && N_old_state[point.idx] == 0
                point.t = time
            end
        end
    end
end

#Funkcja zwracająca wektor odległości od obserwatorów węzła o indeksie idx
function getDistanceFromObservers(N::Network, obs::Vector{Observer}, idx::Int)
    d = Vector{Float64}()
    ds = desopo_pape_shortest_paths(N.graph, idx)
    for point in obs
        if ds.dists[point.idx] > floatmax(Float16)
            push!(d, floatmax(Float16))
        else    
            push!(d, float(ds.dists[point.idx]))
        end    
    end
    return d
end

#Funkcja zwracająca wyniki algorytmu korelacyjnego dla wszystkich węzłów sieci
function getScore(N::Network, obs::Vector{Observer})
    score = Vector{Float64}()
    t = Vector{Float64}()
    for point in obs
        push!(t, float(point.t))
    end
    for i in 1:length(N.network_state)
        d = getDistanceFromObservers(N, obs, i)
        sc::Float64 = cor(t, d)
        if isnan(sc)
            sc = -1.0   
        end
        push!(score, sc)
    end
    return score
end

#Funkcja wyznaczająca precyzję i ranking w symulacji na podstawie wektora wyników
function analizeScore(N::Network, score::Vector{Float64})
    solutions = Vector{Int}()

    for i in 1:length(score)
        if abs(score[i] - maximum(score)) < 0.001
            push!(solutions, i)
        end
    end

    src_score = score[N.source_idx]
    if N.source_idx in solutions
        prec = 1.0 / length(solutions)
    else
        prec = 0.0
    end

    rank = maximum(findall(x -> x == src_score, sort(score, rev=true)))
    return prec, rank

end

#Funkcja przeprowadzająca pojedyńczą symulacje 
function algorithm(N::Network, inf_prob_loc::Float64, gamma::Float64, observer_count::Int) 
    obs = getObservers(N, observer_count) 
    all_obs_infected::Bool = false
    time_step::Int = 1
    
    while all_obs_infected == false
        N_temp_vect = copy(N.network_state)
        N = @set N.network_state = get_next_step(N, inf_prob_loc, gamma)
        actuateObservers(N::Network, N_temp_vect, obs, time_step)
        all_obs_infected = true
        
        for observer in obs
            if observer.t == Int(floatmax(Float16)) 
                all_obs_infected = false
            end
        end    
        time_step += 1
        if time_step > 100000   #Warunek brzegowy symulacji
            println("Utknieto w pętli")
            break
        end
    end
   
    prec_kor, rank_kor = analizeScore(N, getScore(N, obs)) 
    return prec_kor, rank_kor, obs, N
      
end

#Funkcja resetująca stan sieci
function resetExistingNetwork(N::Network)
    new_network_state = Vector{Int}()
    for i in 1:length(N.network_state)
        if i == N.source_idx
            push!(new_network_state, 1)
        else
            push!(new_network_state, 0)
        end
    end
    N = @set N.network_state = new_network_state
end


#Funkcja przeprowadzająca serię symulacji w trybie beta i zapisująca wyniki symulacji do pliku data.txt
function main_beta(gamma::Float64, network_params, observer_count::Int, beta_start::Float64, beta_step::Float64, i_max::Int, j_max::Int)
    file = open("data.txt", "w")
    for i in 1:i_max
        rank_avg_vect_kor = Vector{Float64}()
        prec_avg_vect_kor = Vector{Float16}()
        beta_vect = Vector{Float64}()
        beta = beta_start + beta_step * (i - 1)
        push!(beta_vect, beta)
        for j in 1:j_max

            if length(network_params) == 2
                nodes::Int = network_params[1]
                prob::Float16 = network_params[2]
                N::Network = generate_network(nodes, prob) 
            elseif length(network_params) == 3
                nodes = network_params[1]
                base_nodes::Int = network_params[2]
                edges::Int = network_params[3]
                N = generate_network(nodes, base_nodes, edges) 
            end
            
            prec_kor, rank_kor = algorithm(N, beta, gamma, observer_count)
            push!(prec_avg_vect_kor, prec_kor)
            push!(rank_avg_vect_kor, rank_kor)
            resetExistingNetwork(N)
              
        end
        
        prec_avg_kor = sum(prec_avg_vect_kor) / length(prec_avg_vect_kor)
        std_dev_prec_kor = std(prec_avg_vect_kor) / length(prec_avg_vect_kor)
        rank_avg_kor = sum(rank_avg_vect_kor) / length(rank_avg_vect_kor)
        std_dev_rank_kor = std(rank_avg_vect_kor) / length(rank_avg_vect_kor)
        println("srednia Precyzja (korelacyjny):  ", prec_avg_kor, " +/-: ", std_dev_prec_kor)
        println("sredni Ranking (korelacyjny):  ", rank_avg_kor, " +/-: ", std_dev_rank_kor)

        write(file, string(beta) * " " * string(prec_avg_kor) * " " * string(rank_avg_kor) * " " * string(std_dev_prec_kor) * " " * string(std_dev_rank_kor) * " \n")
           
    end
    close(file)
end

#Funkcja przeprowadzająca serię symulacji w trybie gamma i zapisująca wyniki symulacji do pliku data.txt
function main_gamma(beta::Float64, network_params, observer_count::Int, gamma_start::Float64, gamma_step::Float64, i_max::Int, j_max::Int)
    file = open("data.txt", "w")
    for i in 1:i_max
        rank_avg_vect_kor = Vector{Float64}()
        prec_avg_vect_kor = Vector{Float16}()
        gamma_vect = Vector{Float64}()
        gamma = gamma_start + gamma_step * (i - 1)  
        push!(gamma_vect, gamma)
        for j in 1:j_max

            if length(network_params) == 2
                nodes::Int = network_params[1]
                prob::Float16 = network_params[2]
                N::Network = generate_network(nodes, prob) 
            elseif length(network_params) == 3
                nodes = network_params[1]
                base_nodes::Int = network_params[2]
                edges::Int = network_params[3]
                N = generate_network(nodes, base_nodes, edges) 
            end
                prec_kor, rank_kor = algorithm(N, beta, gamma, observer_count)
                push!(prec_avg_vect_kor, prec_kor)
                push!(rank_avg_vect_kor, rank_kor)
                resetExistingNetwork(N)
        end
            prec_avg_kor = sum(prec_avg_vect_kor) / length(prec_avg_vect_kor)
            std_dev_prec_kor = std(prec_avg_vect_kor) / length(prec_avg_vect_kor)
            rank_avg_kor = sum(rank_avg_vect_kor) / length(rank_avg_vect_kor)
            std_dev_rank_kor = std(rank_avg_vect_kor) / length(rank_avg_vect_kor)
            println("srednia Precyzja (korelacyjny):  ", prec_avg_kor, " +/-: ", std_dev_prec_kor)
            println("sredni Ranking (korelacyjny):  ", rank_avg_kor, " +/-: ", std_dev_rank_kor)

            write(file, string(gamma) * " " * string(prec_avg_kor) * " " * string(rank_avg_kor) * " " * string(std_dev_prec_kor) * " " * string(std_dev_rank_kor) * " \n")
    end
    close(file)
end


#Wczytywanie parametrów serii symulacji z pliku params.txt
params_array = Vector{}()
file = open("params.txt", "r")
for line in readlines(file)
    data = split(split(line, "#")[1], " ")
    push!(params_array, data)
end  
close(file) 
if length(params_array[1]) == 2
    erdos_reyni::Bool = true
    n = parse(Int, params_array[1][1])
    p = parse(Float64, params_array[1][2])
    network_params = [n, p]
elseif length(params_array[1]) == 3
    erdos_reyni = false
    n = parse(Int, params_array[1][1])
    n0 = parse(Int, params_array[1][2])
    k = parse(Int, params_array[1][3])
    network_params = [n, n0, k]
end

observer_count::Int = parse(Int, params_array[2][1])
if params_array[3][1] == "gamma"
    use_gamma::Bool = true
    beta = parse(Float64, params_array[4][1])
    gamma_start = parse(Float64, params_array[5][1])
    gamma_step = parse(Float64, params_array[5][2])
    i_max = parse(Int, params_array[5][3])
elseif params_array[3][1] == "beta"
    use_gamma = false
    beta_start = parse(Float64, params_array[4][1])
    beta_step = parse(Float64, params_array[4][2])
    i_max = parse(Int, params_array[4][3])
    gamma = parse(Float64, params_array[5][1])
end

j_max::Int = parse(Int, params_array[6][1])

#Wywołanie odpowiedniej serii symulacji 
if use_gamma == true
    main_gamma(beta, network_params, observer_count, gamma_start, gamma_step, i_max, j_max)
elseif use_gamma == false
    main_beta(gamma, network_params, observer_count, beta_start, beta_step, i_max, j_max)    
end

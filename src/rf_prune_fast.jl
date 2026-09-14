using LinearAlgebra

if !isdefined(@__MODULE__, :RFSparseCladoEvents)
    include(joinpath(@__DIR__, "rf_prune_cache.jl"))
end

function Apartial_ll_fast!(start_states::Vector{Float64},
                           end_states::Vector{Float64},
                           QPmat::Matrix{Float64})
    mul!(start_states, QPmat, end_states)
    return start_states
end

function Cpartial_ll_sparse!(parent_states::Vector{Float64},
                             left_child_states::Vector{Float64},
                             right_child_states::Vector{Float64},
                             clado_events::RFSparseCladoEvents)
    for parent in eachindex(parent_states)
        total = 0.0
        for event in clado_events.by_parent[parent]
            total += event.prob * left_child_states[event.left] * right_child_states[event.right]
        end
        parent_states[parent] = total
    end

    return parent_states
end

function rf_prune_algo_fast(Q::Matrix{Float64},
                            tip_probs::Vector{Vector{Float64}},
                            cladoPmat::Array{Float64, 3},
                            brs::Matrix{Float64},
                            node_path::Vector{Vector{Int64}})
    clado_events = make_sparse_clado_events(cladoPmat)
    transition_cache = make_branch_transition_cache(Q, brs, node_path)

    return rf_prune_algo_fast_cached(tip_probs, clado_events, transition_cache, node_path)
end

function rf_prune_algo_fast_cached(tip_probs::Vector{Vector{Float64}},
                                   clado_events::RFSparseCladoEvents,
                                   transition_cache::RFBranchTransitionCache,
                                   node_path::Vector{Vector{Int64}})
    uf_thresh = 10^-100
    uf_counter = 0

    ntips = length(tip_probs)
    nstates = clado_events.nstates

    edge_probs = vcat(tip_probs, [fill(0.0, nstates) for _ in 1:(ntips - 1)])
    post_clado_probs = vcat(tip_probs, [fill(0.0, nstates) for _ in 1:(ntips - 1)])

    Ll = zeros(nstates)
    Lr = zeros(nstates)
    parent_states = zeros(nstates)

    for trav_node in node_path
        parent = trav_node[1]
        left_child = trav_node[2]
        right_child = trav_node[3]

        Pt_l = transition_cache.by_child[left_child]
        Pt_r = transition_cache.by_child[right_child]
        Pt_l === nothing && error("No cached transition matrix for child node $(left_child)")
        Pt_r === nothing && error("No cached transition matrix for child node $(right_child)")

        Apartial_ll_fast!(Ll, edge_probs[left_child], Pt_l)
        if any(<(uf_thresh), Ll)
            Ll .*= 10.0^100.0
            uf_counter += 1
        end

        Apartial_ll_fast!(Lr, edge_probs[right_child], Pt_r)
        if any(<(uf_thresh), Lr)
            Lr .*= 10.0^100.0
            uf_counter += 1
        end

        post_clado_probs[right_child] = copy(Lr)
        post_clado_probs[left_child] = copy(Ll)

        Cpartial_ll_sparse!(parent_states, Ll, Lr, clado_events)
        edge_probs[parent] = copy(parent_states)
    end

    logL = log(sum(edge_probs[ntips + 1])) + (-100 * uf_counter * log(10))

    return logL::Float64, edge_probs::Vector{Vector{Float64}}, post_clado_probs::Vector{Vector{Float64}}
end

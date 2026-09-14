if !isdefined(@__MODULE__, :rf_prune_algo_fast_cached)
    include(joinpath(@__DIR__, "rf_prune_fast.jl"))
end

if !isdefined(@__MODULE__, :gillespie_rejection_map)
    include(joinpath(@__DIR__, "rf_stochastic_mapping.jl"))
end

struct RFCladoSupportEvent
    left::Int
    right::Int
end

struct RFMCMCFastPruneCache
    brs::Matrix{Float64}
    node_path::Vector{Vector{Int64}}
    branch_lengths::Vector{Union{Nothing, Float64}}
    clado_support::Vector{Vector{RFCladoSupportEvent}}
    nstates::Int
    nnz::Int
end

function make_sparse_clado_support(cladoPmat::Array{Float64, 3})
    nparents, nleft, nright = size(cladoPmat)
    nleft == nparents || error("cladoPmat left-child dimension does not match parent dimension")
    nright == nparents || error("cladoPmat right-child dimension does not match parent dimension")

    by_parent = [RFCladoSupportEvent[] for _ in 1:nparents]
    nnz = 0

    for parent in 1:nparents
        parent_support = by_parent[parent]
        for left in 1:nleft
            for right in 1:nright
                if cladoPmat[parent, left, right] != 0.0
                    push!(parent_support, RFCladoSupportEvent(left, right))
                    nnz += 1
                end
            end
        end
    end

    return by_parent, nnz
end

function make_rf_mcmc_fast_prune_cache(brs::Matrix{Float64},
                                       node_path::Vector{Vector{Int64}},
                                       cladoPmat::Array{Float64, 3})
    branch_lengths = make_child_branch_lengths(brs, node_path)
    clado_support, nnz = make_sparse_clado_support(cladoPmat)
    return RFMCMCFastPruneCache(brs, node_path, branch_lengths, clado_support, size(cladoPmat, 1), nnz)
end

function make_rf_mcmc_fast_prune_cache(tree::rtree,
                                       cladoPmat::Array{Float64, 3};
                                       rescale_tree::Bool = true)
    ntips = length(tree.tlab)
    edges = cat(tree.ed, [2 * ntips ntips + 1], dims = 1)
    triads = maketriads(edges)
    node_path = get_trav_path_4prun(triads)
    brs = sortslices(branching_times(tree), dims = 1, by = x -> x[5], rev = true)

    if rescale_tree
        brs[:, 3:end] = brs[:, 3:end] ./ maximum(brs[:, 3:end])
    end

    return make_rf_mcmc_fast_prune_cache(brs, node_path, cladoPmat)
end

function make_sparse_clado_events(cache::RFMCMCFastPruneCache,
                                  cladoPmat::Array{Float64, 3})
    size(cladoPmat, 1) == cache.nstates || error("cladoPmat state count does not match cache")
    by_parent = [RFCladoEvent[] for _ in 1:cache.nstates]
    nnz = 0

    for parent in 1:cache.nstates
        parent_events = by_parent[parent]
        for support_event in cache.clado_support[parent]
            prob = cladoPmat[parent, support_event.left, support_event.right]
            if prob != 0.0
                push!(parent_events, RFCladoEvent(support_event.left, support_event.right, prob))
                nnz += 1
            end
        end
    end

    return RFSparseCladoEvents(by_parent, nnz, cache.nstates)
end

function make_branch_transition_cache(Q::Matrix{Float64},
                                      cache::RFMCMCFastPruneCache)
    return make_branch_transition_cache(Q, cache.branch_lengths, cache.node_path)
end

function rf_prune_algo_fast_mcmc(Q::Matrix{Float64},
                                 tip_probs::Vector{Vector{Float64}},
                                 cladoPmat::Array{Float64, 3},
                                 cache::RFMCMCFastPruneCache)
    transition_cache = make_branch_transition_cache(Q, cache)
    clado_events = make_sparse_clado_events(cache, cladoPmat)
    return rf_prune_algo_fast_cached(tip_probs, clado_events, transition_cache, cache.node_path)
end

function rf_prune_algo_fast_mcmc(tip_probs::Vector{Vector{Float64}},
                                 cladoPmat::Array{Float64, 3},
                                 transition_cache::RFBranchTransitionCache,
                                 cache::RFMCMCFastPruneCache)
    clado_events = make_sparse_clado_events(cache, cladoPmat)
    return rf_prune_algo_fast_cached(tip_probs, clado_events, transition_cache, cache.node_path)
end

function rate_par_upd_fast(rate_pars_c,
                           tip_probs,
                           clado_Pmat_c,
                           lL_c,
                           prunelL_c,
                           clado_priorlL_c,
                           rate_priorlL_c,
                           edge_probs,
                           post_clado_probs,
                           Q_par_matrix,
                           Q_zeros,
                           Q_index_vec,
                           prior_dists,
                           transition_cache_c::RFBranchTransitionCache,
                           fast_cache::RFMCMCFastPruneCache,
                           tuning_par_vec,
                           par_acceptfreq_vec,
                           par_propfreq_vec,
                           prior_only = false)
    rate_pars_p = copy(rate_pars_c)
    lL_p = copy(lL_c)
    prunelL_p = copy(prunelL_c)
    rate_priorlL_p = copy(rate_priorlL_c)
    transition_cache_p = transition_cache_c

    prop_par_ind = rand(eachindex(rate_pars_c))

    if rate_pars_c[prop_par_ind] == 0.0000
        lL_ratio = -1.0
    else
        rate_pars_p[prop_par_ind], hr = multi_move(rate_pars_c[prop_par_ind],
                                                   tuning_par_vec[prop_par_ind])

        Q_p = fill_emptyQmat(Q_zeros, rate_pars_p, Q_index_vec, Q_par_matrix)
        transition_cache_p = make_branch_transition_cache(Q_p, fast_cache)

        prunelL_p, edge_probs_p, post_clado_probs_p =
            rf_prune_algo_fast_mcmc(tip_probs, clado_Pmat_c, transition_cache_p, fast_cache)

        rate_priorlL_p = sum([log(pdf(prior_dists[i], rate_pars_p[i])) for i in eachindex(rate_pars_p)])

        if prior_only
            lL_p = rate_priorlL_p + clado_priorlL_c + hr
        else
            lL_p = rate_priorlL_p + clado_priorlL_c + prunelL_p + hr
        end

        lL_ratio = exp(lL_p - lL_c)
    end

    if rand(Uniform(0, 1)) < lL_ratio
        rate_pars_c = copy(rate_pars_p)
        rate_priorlL_c = copy(rate_priorlL_p)
        prunelL_c = copy(prunelL_p)
        lL_c = copy(lL_p)
        edge_probs = edge_probs_p
        post_clado_probs = post_clado_probs_p
        transition_cache_c = transition_cache_p

        par_acceptfreq_vec[prop_par_ind] = par_acceptfreq_vec[prop_par_ind] .+ 1
    end

    par_propfreq_vec[prop_par_ind] = par_propfreq_vec[prop_par_ind] .+ 1

    return (rate_pars_c,
            edge_probs,
            post_clado_probs,
            rate_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            par_acceptfreq_vec,
            par_propfreq_vec)
end

function clado_par_upd_fast(rate_pars_c,
                            clado_probs_c,
                            clado_Pmat_c,
                            lL_c,
                            prunelL_c,
                            clado_priorlL_c,
                            rate_priorlL_c,
                            edge_probs,
                            post_clado_probs,
                            clado_prior,
                            split_vec,
                            sub_vec,
                            Q_par_matrix,
                            Q_zeros,
                            Q_index_vec,
                            tip_probs,
                            transition_cache_c::RFBranchTransitionCache,
                            fast_cache::RFMCMCFastPruneCache,
                            clado_tuning_par_vec,
                            par_acceptfreq_vec,
                            par_propfreq_vec,
                            prior_only = false)
    clado_probs_p = copy(clado_probs_c)
    lL_p = copy(lL_c)
    prunelL_p = copy(prunelL_c)
    clado_priorlL_p = copy(clado_priorlL_c)

    par_up = rand(1:2)
    clado_probs_p[par_up], hr = clado_multi_move(clado_probs_c[par_up], clado_tuning_par_vec[1])
    clado_probs_p[3] = 1 - sum(clado_probs_p[1:2])

    if 0 < sum(clado_probs_p[1:2]) < 1.0
        clado_Pmat_p = fill_subsplit_cladoPmat(clado_Pmat_c,
                                               clado_probs_p,
                                               split_vec,
                                               sub_vec)

        prunelL_p, edge_probs_p, post_clado_probs_p =
            rf_prune_algo_fast_mcmc(tip_probs, clado_Pmat_p, transition_cache_c, fast_cache)

        clado_priorlL_p = sum([log(pdf(clado_prior[i], clado_probs_p[i])) for i in eachindex(clado_prior)])

        if prior_only
            lL_p = clado_priorlL_p + rate_priorlL_c + hr
        else
            lL_p = clado_priorlL_p + rate_priorlL_c + prunelL_p + hr
        end

        lL_ratio = exp(lL_p - lL_c)

        if rand(Uniform(0, 1)) < lL_ratio
            clado_probs_c = copy(clado_probs_p)
            clado_priorlL_c = copy(clado_priorlL_p)
            prunelL_c = copy(prunelL_p)
            lL_c = copy(lL_p)
            edge_probs = edge_probs_p
            post_clado_probs = post_clado_probs_p

            par_acceptfreq_vec[[par_up, 3]] = par_acceptfreq_vec[[par_up, 3]] .+ 1
        end
    end

    par_propfreq_vec[[par_up, 3]] = par_propfreq_vec[[par_up, 3]] .+ 1

    return (clado_probs_c,
            edge_probs,
            post_clado_probs,
            clado_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            par_acceptfreq_vec,
            par_propfreq_vec)
end

function rjswitch_rate_par_upd_fast(rate_pars_c,
                                    rate_pars_switches_c,
                                    tip_probs,
                                    clado_Pmat_c,
                                    lL_c,
                                    prunelL_c,
                                    clado_priorlL_c,
                                    rate_priorlL_c,
                                    edge_probs,
                                    post_clado_probs,
                                    Q_par_matrix,
                                    Q_zeros,
                                    Q_index_vec,
                                    prior_dists,
                                    transition_cache_c::RFBranchTransitionCache,
                                    fast_cache::RFMCMCFastPruneCache,
                                    par_acceptfreq_vec,
                                    par_propfreq_vec,
                                    prior_only = false)
    rate_pars_p = copy(rate_pars_c)
    lL_p = copy(lL_c)
    prunelL_p = copy(prunelL_c)
    rate_priorlL_p = copy(rate_priorlL_c)
    rate_pars_switches_p = copy(rate_pars_switches_c)

    prop_par_ind = rand(eachindex(rate_pars_c))

    rate_pars_p[prop_par_ind], hr = birth_death_move(rate_pars_c[prop_par_ind],
                                                     prior_dists[prop_par_ind])
    rate_pars_switches_p = rate_pars_p .!= 0.0

    Q_p = fill_emptyQmat(Q_zeros, rate_pars_p, Q_index_vec, Q_par_matrix)
    transition_cache_p = make_branch_transition_cache(Q_p, fast_cache)

    prunelL_p, edge_probs_p, post_clado_probs_p =
        rf_prune_algo_fast_mcmc(tip_probs, clado_Pmat_c, transition_cache_p, fast_cache)

    rate_priorlL_p = sum([log(pdf(prior_dists[i], rate_pars_p[i])) for i in eachindex(rate_pars_p)])

    if prior_only
        lL_p = rate_priorlL_p + clado_priorlL_c + hr
    else
        lL_p = rate_priorlL_p + clado_priorlL_c + prunelL_p + hr
    end

    lL_ratio = exp(lL_p - lL_c)

    if rand(Uniform(0, 1)) < lL_ratio
        rate_pars_c = copy(rate_pars_p)
        rate_pars_switches_c = copy(rate_pars_switches_p)
        rate_priorlL_c = copy(rate_priorlL_p)
        prunelL_c = copy(prunelL_p)
        lL_c = copy(lL_p)
        edge_probs = edge_probs_p
        post_clado_probs = post_clado_probs_p
        transition_cache_c = transition_cache_p

        par_acceptfreq_vec[prop_par_ind] = par_acceptfreq_vec[prop_par_ind] .+ 1
    end

    par_propfreq_vec[prop_par_ind] = par_propfreq_vec[prop_par_ind] .+ 1

    return (rate_pars_c,
            rate_pars_switches_c,
            edge_probs,
            post_clado_probs,
            rate_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            par_acceptfreq_vec,
            par_propfreq_vec)
end

function rjswitch_clado_par_upd_fast(rate_pars_c,
                                     clado_probs_c,
                                     clado_probs_switches_c,
                                     clado_Pmat_c,
                                     lL_c,
                                     prunelL_c,
                                     clado_priorlL_c,
                                     rate_priorlL_c,
                                     edge_probs,
                                     post_clado_probs,
                                     clado_prior,
                                     split_vec,
                                     sub_vec,
                                     Q_par_matrix,
                                     Q_zeros,
                                     Q_index_vec,
                                     tip_probs,
                                     transition_cache_c::RFBranchTransitionCache,
                                     fast_cache::RFMCMCFastPruneCache,
                                     par_acceptfreq_vec,
                                     par_propfreq_vec,
                                     prior_only = false)
    clado_probs_p = copy(clado_probs_c)
    lL_p = copy(lL_c)
    prunelL_p = copy(prunelL_c)
    clado_priorlL_p = copy(clado_priorlL_c)

    par_up = rand(1:3)
    clado_probs_p[par_up], hr = birth_death_move(clado_probs_c[par_up], clado_prior[1])

    clado_probs_p = clado_probs_p ./ sum(clado_probs_p)
    clado_probs_switches_p = clado_probs_p .!= 0.0

    clado_Pmat_p = fill_subsplit_cladoPmat(clado_Pmat_c,
                                           clado_probs_p,
                                           split_vec,
                                           sub_vec)

    prunelL_p, edge_probs_p, post_clado_probs_p =
        rf_prune_algo_fast_mcmc(tip_probs, clado_Pmat_p, transition_cache_c, fast_cache)

    clado_priorlL_p = sum([log(pdf(clado_prior[i], clado_probs_p[i])) for i in eachindex(clado_prior)])

    if prior_only
        lL_p = clado_priorlL_p + rate_priorlL_c + hr
    else
        lL_p = clado_priorlL_p + rate_priorlL_c + prunelL_p + hr
    end

    lL_ratio = exp(lL_p - lL_c)

    if rand(Uniform(0, 1)) < lL_ratio
        clado_probs_c = copy(clado_probs_p)
        clado_probs_switches_c = copy(clado_probs_switches_p)
        clado_priorlL_c = copy(clado_priorlL_p)
        prunelL_c = copy(prunelL_p)
        lL_c = copy(lL_p)
        edge_probs = edge_probs_p
        post_clado_probs = post_clado_probs_p

        par_acceptfreq_vec[par_up] = par_acceptfreq_vec[par_up] .+ 1
    end

    par_propfreq_vec[[par_up]] = par_propfreq_vec[[par_up]] .+ 1

    return (clado_probs_c,
            clado_probs_switches_c,
            edge_probs,
            post_clado_probs,
            clado_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            par_acceptfreq_vec,
            par_propfreq_vec)
end

function realfun_mcmc_fast(start_rates, start_clado_probs,
                           iters, iter_trims, write_interval,
                           anc_state_sampling,
                           Q_par_matrix, Q_index_vec, move_types, rate_prior_dists,
                           clado_Pmat_unpar, clado_prior_dists,
                           rf_states,
                           tip_probs,
                           tree,
                           log_filename,
                           rate_tuning_par, clado_tuning_par,
                           proposal_probs,
                           prior_only = false,
                           rescale_tree = true,
                           write_simmap = false,
                           simmap_max_attempts = 10_000,
                           simmap_seed = 90210)
    ed = tree.ed
    ntips = length(tree.tlab)

    edges = cat(ed, [2 * ntips ntips + 1], dims = 1)
    triads = maketriads(edges)
    node_path = get_trav_path_4prun(triads)
    brs = sortslices(branching_times(tree), dims = 1, by = x -> x[5], rev = true)

    if rescale_tree
        brs[:, 3:end] = brs[:, 3:end] / maximum(brs[:, 3:end])
    end

    root_probs = fill(1.0, length(rf_states))
    Q_zeros = zeros(length(rf_states), length(rf_states))
    rate_pars_c = copy(start_rates)
    Q_start = fill_emptyQmat(Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)

    clado_probs_c = copy(start_clado_probs[:, 2])
    _, split_index_vec, sub_index_vec = get_cladoPmat_par_vecs(clado_Pmat_unpar)
    clado_Pmat_support = copy(clado_Pmat_unpar)
    clado_Pmat_c = fill_subsplit_cladoPmat(clado_Pmat_unpar,
                                           start_clado_probs[:, 2],
                                           split_index_vec,
                                           sub_index_vec)

    rate_pars_switches_c = fill(true, length(rate_pars_c))
    clado_probs_switches_c = fill(true, length(clado_probs_c))

    rate_tuning_par_vec = [rate_tuning_par for _ in move_types]
    rate_pars_acceptfreq_vec = [0 for _ in rate_pars_c]
    rate_pars_propfreq_vec = [0 for _ in rate_pars_c]
    rate_zeroswitch_acceptfreq_vec = [0 for _ in rate_pars_c]
    rate_zeroswitch_propfreq_vec = [0 for _ in rate_pars_c]

    clado_tuning_par_vec = [clado_tuning_par for _ in clado_probs_c]
    clado_par_acceptfreq_vec = [0 for _ in clado_probs_c]
    clado_par_propfreq_vec = [0 for _ in clado_probs_c]
    clado_zeroswitch_acceptfreq_vec = [0 for _ in clado_probs_c]
    clado_zeroswitch_propfreq_vec = [0 for _ in clado_probs_c]

    AR_move_types = ["AR" * move for move in move_types]
    AR_switch_move_types = ["AR" * "_switch" * move for move in move_types]

    log_header_string = "iter" * "\t" *
                        "lL_c" * "\t" *
                        "prunelL_c" * "\t" *
                        "rate_priorlL_c" * "\t" *
                        "clado_priorlL_c" * "\t" *
                        join(move_types, "\t") * "\t" *
                        join(string.("switch", move_types), "\t") * "\t" *
                        join(start_clado_probs[:, 1], "\t") * "\t" *
                        join(string.("switch_", start_clado_probs[:, 1]), "\t") * "\t" *
                        join(AR_move_types, "\t") * "\t" *
                        join(AR_switch_move_types, "\t") * "\t" *
                        "AR_clado_split_prob" * "\t" *
                        "AR_clado_sub_prob" * "\t" *
                        "AR_clado_equal_prob" * "\t" *
                        "AR_switch_clado_split_prob" * "\t" *
                        "AR_switch_clado_sub_prob" * "\t" *
                        "AR_switch_clado_equal_prob"

    anc_states_log_filename = log_filename * "_anc_states"
    anc_clado_log_filename = log_filename * "_anc_clados"
    simmap_log_filename = log_filename * "_stoch.log"

    chain_cache = []
    anc_clado_cache = []
    anc_states_cache = []
    simmap_cache = []

    log_file = open(log_filename, "a")
    println(log_file, log_header_string)
    print(log_header_string)
    close(log_file)

    simmap_buffer = write_simmap ? build_simmap_buffer(tree) : SimmapBranch[]
    simmap_edges = write_simmap ? build_simmap_edges(tree, brs) : SimmapEdge[]
    simmap_rng = Random.MersenneTwister(simmap_seed)

    if write_simmap
        simmap_log_file = open(simmap_log_filename, "a")
        println(simmap_log_file, simmap_log_header(tree))
        close(simmap_log_file)
    end

    fast_cache = make_rf_mcmc_fast_prune_cache(brs, node_path, clado_Pmat_support)
    transition_cache_c = make_branch_transition_cache(Q_start, fast_cache)
    prunelL_c, edge_probs, post_clado_probs =
        rf_prune_algo_fast_mcmc(tip_probs, clado_Pmat_c, transition_cache_c, fast_cache)

    rate_priorlL_c = sum([log(pdf(rate_prior_dists[i], rate_pars_c[i])) for i in eachindex(rate_pars_c)])
    clado_priorlL_c = sum([log(pdf(clado_prior_dists[i], clado_probs_c[i])) for i in eachindex(clado_prior_dists)])

    if prior_only
        lL_c = rate_priorlL_c + clado_priorlL_c
    else
        lL_c = rate_priorlL_c + clado_priorlL_c + prunelL_c
    end

    log_line_string = "0" * "\t" *
                      string(round(lL_c, digits = 4)) * "\t" *
                      string(round(prunelL_c, digits = 4)) * "\t" *
                      string(round(rate_priorlL_c, digits = 4)) * "\t" *
                      string(round(clado_priorlL_c, digits = 4)) * "\t" *
                      join(string.(round.(rate_pars_c, digits = 5)), "\t") * "\t" *
                      join(string.(round.(rate_pars_switches_c, digits = 5)), "\t") * "\t" *
                      join(string.(round.(clado_probs_c, digits = 5)), "\t") * "\t" *
                      join(string.(round.(clado_probs_switches_c, digits = 5)), "\t") * "\t" *
                      join(string.(round.(rate_pars_acceptfreq_vec ./ rate_pars_propfreq_vec, digits = 3)), "\t") * "\t" *
                      join(string.(round.(rate_zeroswitch_acceptfreq_vec ./ rate_zeroswitch_propfreq_vec, digits = 3)), "\t") * "\t" *
                      join(string.(round.(clado_par_acceptfreq_vec ./ clado_par_propfreq_vec, digits = 3)), "\t") * "\t" *
                      join(string.(round.(clado_zeroswitch_acceptfreq_vec ./ clado_zeroswitch_propfreq_vec, digits = 3)), "\t")

    log_file = open(log_filename, "a")
    println(log_file, log_line_string)
    print(log_line_string)
    close(log_file)

    move_names = collect(keys(proposal_probs))
    move_probabilities = collect(values(proposal_probs) ./ sum(values(proposal_probs)))

    for iter in 1:iters
        par2upd = sample(move_names, Weights(move_probabilities))

        if par2upd == :rate_move
            rate_pars_c,
            edge_probs,
            post_clado_probs,
            rate_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            rate_pars_acceptfreq_vec,
            rate_pars_propfreq_vec = rate_par_upd_fast(rate_pars_c,
                                                       tip_probs,
                                                       clado_Pmat_c,
                                                       lL_c,
                                                       prunelL_c,
                                                       clado_priorlL_c,
                                                       rate_priorlL_c,
                                                       edge_probs,
                                                       post_clado_probs,
                                                       Q_par_matrix,
                                                       Q_zeros,
                                                       Q_index_vec,
                                                       rate_prior_dists,
                                                       transition_cache_c,
                                                       fast_cache,
                                                       rate_tuning_par_vec,
                                                       rate_pars_acceptfreq_vec,
                                                       rate_pars_propfreq_vec,
                                                       prior_only)
        elseif par2upd == :clado_move
            clado_probs_c,
            edge_probs,
            post_clado_probs,
            clado_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            clado_par_acceptfreq_vec,
            clado_par_propfreq_vec = clado_par_upd_fast(rate_pars_c,
                                                        clado_probs_c,
                                                        clado_Pmat_c,
                                                        lL_c,
                                                        prunelL_c,
                                                        clado_priorlL_c,
                                                        rate_priorlL_c,
                                                        edge_probs,
                                                        post_clado_probs,
                                                        clado_prior_dists,
                                                        split_index_vec,
                                                        sub_index_vec,
                                                        Q_par_matrix,
                                                        Q_zeros,
                                                        Q_index_vec,
                                                        tip_probs,
                                                        transition_cache_c,
                                                        fast_cache,
                                                        clado_tuning_par_vec,
                                                        clado_par_acceptfreq_vec,
                                                        clado_par_propfreq_vec,
                                                        prior_only)
        elseif par2upd == :rate_zero_switch
            rate_pars_c,
            rate_pars_switches_c,
            edge_probs,
            post_clado_probs,
            rate_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            rate_zeroswitch_acceptfreq_vec,
            rate_zeroswitch_propfreq_vec = rjswitch_rate_par_upd_fast(rate_pars_c,
                                                                      rate_pars_switches_c,
                                                                      tip_probs,
                                                                      clado_Pmat_c,
                                                                      lL_c,
                                                                      prunelL_c,
                                                                      clado_priorlL_c,
                                                                      rate_priorlL_c,
                                                                      edge_probs,
                                                                      post_clado_probs,
                                                                      Q_par_matrix,
                                                                      Q_zeros,
                                                                      Q_index_vec,
                                                                      rate_prior_dists,
                                                                      transition_cache_c,
                                                                      fast_cache,
                                                                      rate_zeroswitch_acceptfreq_vec,
                                                                      rate_zeroswitch_propfreq_vec,
                                                                      prior_only)
        elseif par2upd == :clado_zero_switch
            clado_probs_c,
            clado_probs_switches_c,
            edge_probs,
            post_clado_probs,
            clado_priorlL_c,
            prunelL_c,
            lL_c,
            transition_cache_c,
            clado_zeroswitch_acceptfreq_vec,
            clado_zeroswitch_propfreq_vec = rjswitch_clado_par_upd_fast(rate_pars_c,
                                                                        clado_probs_c,
                                                                        clado_probs_switches_c,
                                                                        clado_Pmat_c,
                                                                        lL_c,
                                                                        prunelL_c,
                                                                        clado_priorlL_c,
                                                                        rate_priorlL_c,
                                                                        edge_probs,
                                                                        post_clado_probs,
                                                                        clado_prior_dists,
                                                                        split_index_vec,
                                                                        sub_index_vec,
                                                                        Q_par_matrix,
                                                                        Q_zeros,
                                                                        Q_index_vec,
                                                                        tip_probs,
                                                                        transition_cache_c,
                                                                        fast_cache,
                                                                        clado_zeroswitch_acceptfreq_vec,
                                                                        clado_zeroswitch_propfreq_vec,
                                                                        prior_only)
        end

        if mod(iter, anc_state_sampling) == 0
            anc_states, anc_clado_events = sample_anc_states(edge_probs,
                                                             post_clado_probs,
                                                             root_probs,
                                                             tip_probs,
                                                             Q_zeros,
                                                             rate_pars_c,
                                                             Q_index_vec,
                                                             Q_par_matrix,
                                                             clado_Pmat_c,
                                                             brs,
                                                             node_path)

            push!(anc_states_cache, join(string.(anc_states), "\t"))
            push!(anc_clado_cache, join(string.(anc_clado_events), "\t"))

            if write_simmap
                Q_simmap = fill_emptyQmat(Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)
                fill_simmap_buffer!(simmap_buffer,
                                     Q_simmap,
                                     anc_states,
                                     simmap_edges;
                                     max_attempts = simmap_max_attempts,
                                     rng = simmap_rng)
                push!(simmap_cache, simmap_log_row(iter, tree, simmap_buffer, simmap_edges))
            end

            if mod(length(anc_states_cache), write_interval) == 0
                anc_states_log_file = open(anc_states_log_filename, "a")
                [println(anc_states_log_file, anc_states_cache[i]) for i in eachindex(anc_states_cache)]
                close(anc_states_log_file)
                anc_states_cache = []
            end

            if mod(length(anc_clado_cache), write_interval) == 0
                anc_clado_log_file = open(anc_clado_log_filename, "a")
                [println(anc_clado_log_file, anc_clado_cache[i]) for i in eachindex(anc_clado_cache)]
                close(anc_clado_log_file)
                anc_clado_cache = []
            end

            if write_simmap && mod(length(simmap_cache), write_interval) == 0
                simmap_log_file = open(simmap_log_filename, "a")
                [println(simmap_log_file, simmap_cache[i]) for i in eachindex(simmap_cache)]
                close(simmap_log_file)
                simmap_cache = []
            end
        end

        if mod(iter, iter_trims) == 0
            file_line = string(iter) * "\t" *
                        string(round(lL_c, digits = 4)) * "\t" *
                        string(round(prunelL_c, digits = 4)) * "\t" *
                        string(round(rate_priorlL_c, digits = 4)) * "\t" *
                        string(round(clado_priorlL_c, digits = 4)) * "\t" *
                        join(string.(round.(rate_pars_c, digits = 5)), "\t") * "\t" *
                        join(string.(round.(rate_pars_switches_c, digits = 5)), "\t") * "\t" *
                        join(string.(round.(clado_probs_c, digits = 5)), "\t") * "\t" *
                        join(string.(round.(clado_probs_switches_c, digits = 5)), "\t") * "\t" *
                        join(string.(round.(rate_pars_acceptfreq_vec ./ rate_pars_propfreq_vec, digits = 3)), "\t") * "\t" *
                        join(string.(round.(rate_zeroswitch_acceptfreq_vec ./ rate_zeroswitch_propfreq_vec, digits = 3)), "\t") * "\t" *
                        join(string.(round.(clado_par_acceptfreq_vec ./ clado_par_propfreq_vec, digits = 3)), "\t") * "\t" *
                        join(string.(round.(clado_zeroswitch_acceptfreq_vec ./ clado_zeroswitch_propfreq_vec, digits = 3)), "\t")

            push!(chain_cache, file_line)

            if mod(length(chain_cache), write_interval) == 0
                log_file = open(log_filename, "a")
                [println(log_file, chain_cache[i]) for i in eachindex(chain_cache)]
                close(log_file)
                chain_cache = []
            end
        end
    end

    if write_simmap && !isempty(simmap_cache)
        simmap_log_file = open(simmap_log_filename, "a")
        [println(simmap_log_file, simmap_cache[i]) for i in eachindex(simmap_cache)]
        close(simmap_log_file)
    end
end

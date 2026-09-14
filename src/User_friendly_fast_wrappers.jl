if !isdefined(@__MODULE__, :realfun_mcmc_fast)
    include(joinpath(@__DIR__, "rf_mcmc_fast.jl"))
end

function RFBS_fast(;data_fn,
                   tree,
                   max_range,
                   adj_matrix,
                   iters = 1000000,
                   iter_trims = iters / 10000,
                   anc_state_sampling = iters / 10000,
                   write_interval = 100,
                   proposal_probs = (rate_move = 10,
                                     clado_move = 5,
                                     rate_zero_switch = 2.0,
                                     clado_zero_switch = 0.0),
                   rate_tuning_par = 1.5,
                   clado_tuning_par = 0.4,
                   exponential_prior_rate = 0.5,
                   allow_double_gains = true,
                   allow_double_losses = false,
                   allow_single_gains = true,
                   allow_single_losses = true,
                   allow_real_switches = true,
                   allow_realfun_switches = false,
                   ecological = true,
                   allopatric = true,
                   by_biome = false,
                   by_rf = true,
                   by_gainloss = true,
                   by_doublesingle = true,
                   write_simmap = false,
                   simmap_max_attempts = 10_000,
                   simmap_seed = 90210,
                   output_dirname = "RFBS_fast_out")
    tip_probs_new, nbiomes, rf_states = read_data_csv(data_fn; adj_matrix = adj_matrix, max_real_biome = max_range)

    rate_prior_vec = [Exponential(exponential_prior_rate)]

    par_matrix, rf_states, move_types, index_vec, _ = make_rf_par_matrix(nbiomes, max_range,
                                                                         allow_real_switches,
                                                                         allow_realfun_switches,
                                                                         allow_double_gains,
                                                                         allow_double_losses,
                                                                         allow_single_gains,
                                                                         allow_single_losses,
                                                                         by_biome,
                                                                         by_rf,
                                                                         by_gainloss,
                                                                         by_doublesingle,
                                                                         false)

    rf_states_mat = reduce(vcat, permutedims.(rf_states))

    cladoPmat_unpar = makeclado_Pmat_equal_prob(rf_states,
                                                ecological,
                                                allopatric)

    _, clado_probs_start, _, _, _, clado_prior_dists = make_clado_Pmat(rf_states,
                                                                       ecological,
                                                                       allopatric)

    rate_prior_dists = make_Prior(rate_prior_vec, move_types)

    start_rates = [rand(rate_prior_dists[i]) for i in eachindex(move_types)]

    if isdir(output_dirname) == false
        mkdir(output_dirname)
    end

    open(joinpath(output_dirname, "state_space"), "w") do io
        writedlm(io, rf_states_mat)
    end

    log_filename = joinpath(output_dirname, "RFBS_log")

    realfun_mcmc_fast(start_rates,
                      clado_probs_start,
                      iters,
                      iter_trims,
                      write_interval,
                      anc_state_sampling,
                      par_matrix,
                      index_vec,
                      move_types,
                      rate_prior_dists,
                      cladoPmat_unpar,
                      clado_prior_dists,
                      rf_states,
                      tip_probs_new,
                      tree,
                      log_filename,
                      rate_tuning_par,
                      clado_tuning_par,
                      proposal_probs,
                      false,
                      false,
                      write_simmap,
                      simmap_max_attempts,
                      simmap_seed)
end

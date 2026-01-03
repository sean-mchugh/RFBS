#using BenchmarkTools
#using KernelDensity

"""
    sliding_window_proposal(x, delta)

Propose a new value for `x` by adding a uniform random deviate
in [-delta, delta].
"""
function sliding_window_proposal(x, delta)
    # Propose x' = x + U(-delta, delta)
    return x + rand(Uniform(-delta, delta))
end

using Random, Distributions


function birth_death_move(θ::Real, 
                          q_b::Distribution;
                          rng=Random.GLOBAL_RNG)
    if θ == 0.0
        # BIRTH: sample a new parameter from q_b
        β = rand(rng, q_b)
        # forward = p_dim * pdf(q_b, β)
        # reverse = p_dim * 1
        # prop_ratio = reverse / forward = 1 / pdf(q_b, β)
        prop_ratio = 1 / pdf(q_b, β)
        return (β, 1, log(prop_ratio))
    else
        # DEATH: remove the existing parameter
        # forward = p_dim * 1
        # reverse = p_dim * pdf(q_b, θ)
        # prop_ratio = reverse / forward = pdf(q_b, θ)
        prop_ratio = pdf(q_b, θ)
        return (0.0, 0, log(prop_ratio))
    end
end




function   rate_par_upd(rate_pars_c      ,
                        tip_probs        ,
                        clado_Pmat_c     ,
                        lL_c             ,
                        prunelL_c        ,
                        clado_priorlL_c  ,
                        rate_priorlL_c   , 
                        edge_probs  ,
                        post_clado_probs,                      
                        Q_par_matrix     ,
                        Q_zeros          ,
                        Q_index_vec  ,
                        prior_dists      ,
                        brs              ,
                        node_path        ,
                        tuning_par_vec   ,
                        par_acceptfreq_vec   ,
                        par_propfreq_vec  ,
                        prior_only=false
                        )

    rate_pars_p=copy(rate_pars_c)
    lL_p       =copy(lL_c       )
    prunelL_p  =copy(prunelL_c  )
    rate_priorlL_p  =copy(rate_priorlL_c  )


    prop_par_ind=rand(eachindex(rate_pars_c))
    
    if rate_pars_c[prop_par_ind] == 0.0000

        lL_ratio = -1.0

    else
        rate_pars_p[prop_par_ind], hr   = multi_move(rate_pars_c[prop_par_ind], #parameter
                                                     tuning_par_vec[prop_par_ind]  #tuning par
                                                )

    


         Q_p= fill_emptyQmat( Q_zeros, rate_pars_p, Q_index_vec, Q_par_matrix)

         prunelL_p, edge_probs_p, post_clado_probs_p = rf_prune_algo(Q_p     ,
                                                         tip_probs        ,
                                                         clado_Pmat_c       ,
                                                         brs              ,
                                                         node_path        
                                                                  )
        
        
        rate_priorlL_p = sum([log(pdf(prior_dists[i], rate_pars_p[i])) for i in eachindex(rate_pars_p)]) 


        if prior_only 

            lL_p= rate_priorlL_p + clado_priorlL_c +hr

        else

            lL_p= rate_priorlL_p + clado_priorlL_c + prunelL_p  +hr

        end 


        lL_ratio =  exp(lL_p-lL_c)
    end

    if rand(Uniform(0, 1)) < lL_ratio 

        #Q_c         =copy(Q_p)
        rate_pars_c  =copy(rate_pars_p)
        rate_priorlL_c    =copy(rate_priorlL_p )
        prunelL_c  =copy(prunelL_p  )
        lL_c         =copy(lL_p )
        edge_probs      = edge_probs_p         
        post_clado_probs= post_clado_probs_p 




        par_acceptfreq_vec[prop_par_ind] = par_acceptfreq_vec[prop_par_ind] .+1
    end 

    par_propfreq_vec[prop_par_ind] = par_propfreq_vec[prop_par_ind] .+1

    return (
            rate_pars_c ,
            edge_probs  ,
            post_clado_probs,
            rate_priorlL_c   ,
            prunelL_c   , 
            lL_c        ,
            par_acceptfreq_vec,
            par_propfreq_vec  )

end


function   rjswitch_rate_par_upd(rate_pars_c ,
                            rate_pars_switches_c,
                            tip_probs        ,
                            clado_Pmat_c     ,
                            lL_c             ,
                            prunelL_c        ,
                            clado_priorlL_c  ,
                            rate_priorlL_c   , 
                            edge_probs       ,
                            post_clado_probs ,                      
                            Q_par_matrix     ,
                            Q_zeros          ,
                            Q_index_vec      ,
                            prior_dists      ,
                            brs              ,
                            node_path        ,
                            par_acceptfreq_vec   ,
                            par_propfreq_vec  ,
                            prior_only=false
                            )

    rate_pars_p=copy(rate_pars_c)
    lL_p       =copy(lL_c       )
    prunelL_p  =copy(prunelL_c  )
    rate_priorlL_p  =copy(rate_priorlL_c  )
    rate_pars_switches_p = copy( rate_pars_switches_c    )


    prop_par_ind=rand(eachindex(rate_pars_c))

    rate_pars_p[prop_par_ind], hr   = birth_death_move(rate_pars_c[prop_par_ind], #parameter
                                                       prior_dists[prop_par_ind] # prior dist
                                            )


    rate_pars_switches_p =  rate_pars_p.!=0.0


     Q_p= fill_emptyQmat( Q_zeros, rate_pars_p, Q_index_vec, Q_par_matrix)

     prunelL_p, edge_probs_p, post_clado_probs_p = rf_prune_algo(Q_p     ,
                                                     tip_probs        ,
                                                     clado_Pmat_c       ,
                                                     brs              ,
                                                     node_path        
                                                              )
    
    
    rate_priorlL_p = sum([log(pdf(prior_dists[i], rate_pars_p[i])) for i in eachindex(rate_pars_p)]) 


    if prior_only 

        lL_p= rate_priorlL_p + clado_priorlL_c +hr

    else

        lL_p= rate_priorlL_p + clado_priorlL_c + prunelL_p  +hr

    end 


    lL_ratio =  exp(lL_p-lL_c)

    if rand(Uniform(0, 1)) < lL_ratio 

        #Q_c         =copy(Q_p)
        rate_pars_c     =copy(rate_pars_p)
        rate_pars_switches_c = copy(rate_pars_switches_p)
        rate_priorlL_c  =copy(rate_priorlL_p )
        prunelL_c       =copy(prunelL_p  )
        lL_c            =copy(lL_p )
        edge_probs      = edge_probs_p         
        post_clado_probs= post_clado_probs_p 

        par_acceptfreq_vec[prop_par_ind] = par_acceptfreq_vec[prop_par_ind] .+1
    end 

    par_propfreq_vec[prop_par_ind] = par_propfreq_vec[prop_par_ind] .+1

    return (
            rate_pars_c ,
            rate_pars_switches_c,
            edge_probs  ,
            post_clado_probs,
            rate_priorlL_c   ,
            prunelL_c   , 
            lL_c        ,
            par_acceptfreq_vec,
            par_propfreq_vec  )

end





function   rjswitch_clado_par_upd(rate_pars_c      ,
                         clado_probs_c     , 
                         clado_probs_switches_c , 
                         clado_Pmat_c           ,
                         lL_c                   ,
                         prunelL_c              ,
                         clado_priorlL_c        ,
                         rate_priorlL_c         ,  
                         edge_probs         ,
                         post_clado_probs   ,
                         clado_prior            , 
                         split_vec              ,
                         sub_vec                ,
                         Q_par_matrix     ,
                         Q_zeros          ,
                         Q_index_vec  ,
                         tip_probs              ,
                         brs                    ,
                         node_path              ,
                         par_acceptfreq_vec     ,
                         par_propfreq_vec       ,
                         prior_only=false
                         )

    clado_probs_p     =copy(clado_probs_c )
    lL_p             =copy(lL_c       )
    prunelL_p        =copy(prunelL_c  )
    clado_priorlL_p  =copy(clado_priorlL_c  )
    #rate_priorlL_p  =copy(rate_priorlL_c  )                      


    #prop_par_ind=rand(eachindex(rate_pars_c))

    #rate_pars_p[prop_par_ind], hr   =     multi_move( rate_pars_c[prop_par_ind], #parameter
    #                                                    tuning_par_vec[prop_par_ind]  #tuning par
    #                                                  )
    #
     #print(Q)

     #choose which clado event prob to propose on (equal will always update with either split or sub)

     par_up=rand(1:3)

     #print(clado_probs_c[par_up])
     clado_probs_p[par_up] , hr   = birth_death_move(clado_probs_c[par_up], #parameter
                                               clado_prior[1] # prior dist
                                                    )   


    clado_probs_p = clado_probs_p ./ sum(clado_probs_p)

    clado_probs_switches_p = clado_probs_p .!= 0.0

    print(clado_probs_p)
     #clado_probs_p[3]=1-sum(clado_probs_p[1:2])


    Q= fill_emptyQmat( Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)


     clado_Pmat_p            = fill_subsplit_cladoPmat(clado_Pmat_c, 
                                                       clado_probs_p,
                                                       split_vec,
                                                       sub_vec)

     prunelL_p, edge_probs_p, post_clado_probs_p = rf_prune_algo(Q     ,
                                                            tip_probs        ,
                                                            clado_Pmat_p       ,
                                                            brs              ,
                                                            node_path        
                                                                     )
    

    
    clado_priorlL_p = sum([log(pdf(clado_prior[i], clado_probs_p[i])) for i in eachindex(clado_prior)]) 


    lL_p= clado_priorlL_p + rate_priorlL_c + prunelL_p + hr

    if prior_only 

        lL_p=clado_priorlL_p + rate_priorlL_c + hr

    else

        lL_p=clado_priorlL_p + rate_priorlL_c + prunelL_p + hr

    end 


    lL_ratio =  exp(lL_p - lL_c)

    if rand(Uniform(0, 1)) < lL_ratio 

        clado_probs_c     =copy(clado_probs_p   )
        clado_probs_switches_c     =copy(clado_probs_switches_p   )

        clado_priorlL_c  =copy(clado_priorlL_p)
        prunelL_c        =copy(prunelL_p      )
        lL_c             =copy(lL_p           )
        edge_probs      = edge_probs_p         
        post_clado_probs= post_clado_probs_p 

        par_acceptfreq_vec[par_up] = par_acceptfreq_vec[par_up] .+1
    
    
    end 


    par_propfreq_vec[[par_up]] = par_propfreq_vec[[par_up]] .+1

    return (clado_probs_c       ,
            clado_probs_switches_c       ,
            edge_probs         ,
            post_clado_probs   ,
            clado_priorlL_c          ,
            prunelL_c    , 
            lL_c               ,
            par_acceptfreq_vec ,
            par_propfreq_vec   )

end


#prior_only=false
#rescale_tree=true
#clado_Pmat_unpar = cladoPmat_unpar
#proposal_probs = (rate_move = 10 , clado_move = 5 , rate_zero_switch = 2, clado_zero_switch = 0 )

function realfun_mcmc(start_rates, start_clado_probs, 
                      iters, iter_trims, write_interval,
                      anc_state_sampling, 
                      Q_par_matrix, Q_index_vec, move_types, rate_prior_dists, 
                      clado_Pmat_unpar,  clado_prior_dists,
                      rf_states,
                      tip_probs, 
                      tree, 
                      log_filename, 
                      rate_tuning_par, clado_tuning_par,
                      proposal_probs,
                      prior_only=false ,
                      rescale_tree=true)



    #make prior distribution vectors for rates a cladogenetic splitting probability 
    #rate_prior_dists=make_Prior(rate_prior_vec, move_types)
    #clado_prior_dists=make_Prior(clado_prior_vec, clado_types)
    #print(tree)

    #tip_probs=get_tip_probs(rf_states,tip_states)

    ed=tree.ed
    ntips=length(tree.tlab)
    #writedlm("ets.txt",ets)
    #writedlm("bts.txt", brs)

   # clado_Pmat=makeclado_Pmat(rf_states,rf2stDict)


    edges = cat(ed, [2*ntips ntips + 1], dims = 1)
    triads=maketriads(edges)
    #tip_biomes=tip_areas
    node_path  = get_trav_path_4prun(triads )
    br = branching_times(tree)
    # sort according to branching times
    brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
    
    if(rescale_tree)
        brs[:,3:end]=brs[:,3:end]/maximum(brs[:,3:end])
    end
    root_probs= fill(1.0, length(rf_states))
    #lower = optim_lower_bound
    #upper = optim_upper_bound
    #initial_x = [0.2]
    #inner_optimizer = GradientDescent()
    #start_rates =[rand(prior_dists[i]) for i in eachindex(move_types)]
    #start_rates =[1.0 for i in eachindex(move_types)]



    #generate current Q matrix (doing it at mcmc level instead of ll fn level so it doesnt need to be reassembled every time a par is updated (say a cladogenetic par))
    Q_zeros =zeros(length(rf_states), length(rf_states))
    rate_pars_c=copy(start_rates)
    Q_start= fill_emptyQmat( Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)

    ##generate cladogenetic probability matrix and set as current matrix (either (1) use an unparameterized matrix where all daughter scenarios are equal for any one given parental state, or (2) one where it is weighted by range splitting vs range subsetting

    #clado_split_prob_c=copy(start_clado_split_prob)
    #clado_sub_prob_c=copy(start_clado_sub_prob)
    clado_probs_c=copy(start_clado_probs[:,2])
    String_Clado_Mats, split_index_vec, sub_index_vec=get_cladoPmat_par_vecs(clado_Pmat_unpar)
    #print(clado_probs_c)
    clado_Pmat_c=fill_subsplit_cladoPmat(clado_Pmat_unpar, 
                                         start_clado_probs[:,2],
                                         split_index_vec,
                                         sub_index_vec)

    rate_pars_switches_c = fill(true, length(rate_pars_c))
    clado_probs_switches_c = fill(true, length(clado_probs_c))
                                     

    ###generate tuning parameterrs, aceptance frequencies and proposal frequencies for parameters
    rate_tuning_par_vec      = [rate_tuning_par for i in move_types]
    rate_pars_acceptfreq_vec  = [0 for i in rate_pars_c]
    rate_pars_propfreq_vec    = [0 for i in rate_pars_c]
    rate_zeroswitch_acceptfreq_vec = [0 for i in rate_pars_c]
    rate_zeroswitch_propfreq_vec   = [0 for i in rate_pars_c]


    clado_tuning_par_vec                 =  [clado_tuning_par for i in clado_probs_c]
    clado_par_acceptfreq_vec             =  [0 for i in clado_probs_c]
    clado_par_propfreq_vec               =  [0 for i in clado_probs_c]
    clado_zeroswitch_acceptfreq_vec      =  [0 for i in clado_probs_c]
    clado_zeroswitch_propfreq_vec        =  [0 for i in clado_probs_c]



    AR_move_types=["AR"*move for move in move_types]
    AR_switch_move_types=["AR"* "_switch" * move for move in move_types]

    #log_header_string= "iter" * "\t" * "lL_c" * "\t" * "prunelL_c" * "\t" * "priorlL_c" * "\t" *join( move_types, "\t") * "\t"  * join(  AR_move_types, "\t")
    log_header_string= "iter" * "\t" * 
                        "lL_c" * "\t" *
                        "prunelL_c" * "\t" * 
                        "rate_priorlL_c" * "\t" * 
                        "clado_priorlL_c" * "\t" * 
                        join( move_types, "\t") * "\t"  * 
                        join( string.("switch", move_types), "\t") * "\t"  * 

                        join( start_clado_probs[:,1], "\t") * "\t"  * 
                        join( string.("switch_",start_clado_probs[:,1]), "\t") * "\t"  * 

                        join(  AR_move_types, "\t") *"\t" *
                        join( AR_switch_move_types, "\t") *"\t" *
                        "AR_clado_split_prob" *"\t" *
                        "AR_clado_sub_prob" *"\t" *
                        "AR_clado_equal_prob" * "\t" *
                        "AR_switch_clado_split_prob" *"\t" *
                        "AR_switch_clado_sub_prob" *"\t" *
                        "AR_switch_clado_equal_prob" 





    iter=0

    file_line=log_header_string

    anc_states_log_filename=log_filename*"_anc_states"
    anc_clado_log_filename=log_filename*"_anc_clados"

    #log_filename = open_new_log(log_filename)
    #anc_states_log_filename = open_new_log( anc_st_log_filename)
    #anc_clado_log_filename = open_new_log( anc_clado_log_filename)


    chain_cache=[]
    anc_clado_cache=[]
    anc_states_cache=[]

    #open(log_filename, "w") do log_file

    #    open(anc_states_log_filename, "w") do anc_state_log_file

    #        open(anc_clado_log_filename, "w") do anc_clado_log_file



    log_file =  open(log_filename,"a")
    #write(exampleFileIOStream, file_line)
    println(log_file, file_line)
    print(file_line)

    close(log_file)

    prunelL_c, edge_probs,  post_clado_probs = rf_prune_algo(Q_start      ,
                                                             tip_probs    ,
                                                             clado_Pmat_c ,
                                                             brs          ,
                                                             node_path        
                                                                     )


    rate_priorlL_c = sum([log(pdf(rate_prior_dists[i], rate_pars_c[i])) for i in eachindex(rate_pars_c)]) 

    clado_priorlL_c = sum([log(pdf(clado_prior_dists[i],clado_probs_c[i])) for i in eachindex(clado_prior_dists)]) 

    if prior_only 

        lL_c= rate_priorlL_c + clado_priorlL_c

    else

        lL_c= rate_priorlL_c + clado_priorlL_c + prunelL_c  

    end 



    log_line_string= "0" * "\t" *
    string(round(lL_c, digits=4)) * "\t" * 
    string(round(prunelL_c,digits=4)) * "\t" * 
    string(round(rate_priorlL_c,digits=4)) * "\t" * 
    string(round(clado_priorlL_c, digits=4)) * "\t" *
    join( string.(round.(rate_pars_c, digits=5)), "\t") * "\t"  * 
    join( string.(round.(rate_pars_switches_c, digits=5)), "\t") * "\t"  * 

    join( string.(round.(clado_probs_c, digits=5)), "\t") * "\t"  * 
    join( string.(round.(clado_probs_switches_c, digits=5)), "\t") * "\t"  * 

    #string(round(clado_probs_c, digits=3)) * "\t"  * 
    join( string.(round.(rate_pars_acceptfreq_vec ./rate_pars_propfreq_vec, digits=3)), "\t") * "\t" *
    join( string.(round.(rate_zeroswitch_acceptfreq_vec ./rate_zeroswitch_propfreq_vec, digits=3)), "\t") * "\t" *

    join( string.(round.(clado_par_acceptfreq_vec ./clado_par_propfreq_vec, digits=3)), "\t") * "\t" *
    join( string.(round.(clado_zeroswitch_acceptfreq_vec ./clado_zeroswitch_propfreq_vec, digits=3)), "\t") 


    file_line=log_line_string

    log_file =  open(log_filename,"a")

    println(log_file, file_line)
    print(file_line)

    close(log_file)

   

    #setup proposal prob for proposing on a rate par of clado par

    move_names = collect(keys(proposal_probs))
    move_probabilities = collect(values(proposal_probs)./sum(values(proposal_probs)))



    for iter in 1:iters

        #print("\n")
        #print(iter)
        #print("\n")
       #iter = iter +1
        # Sample from the names using the weights
        par2upd = sample(move_names, Weights(move_probabilities))
        

        #make the first iteration update rate par so Q_c gets made in loop


        if par2upd == :rate_move

            rate_pars_c ,
            edge_probs  , 
            post_clado_probs,
            rate_priorlL_c   ,
            prunelL_c   , 
            lL_c        ,
            rate_pars_acceptfreq_vec ,
            rate_pars_propfreq_vec =rate_par_upd(rate_pars_c      ,
                                                tip_probs        ,
                                                clado_Pmat_c       ,
                                                lL_c             ,
                                                prunelL_c        ,
                                                clado_priorlL_c  ,
                                                rate_priorlL_c   ,  
                                                edge_probs         ,
                                                post_clado_probs   ,                                                
                                                Q_par_matrix     ,
                                                Q_zeros          ,
                                                Q_index_vec  ,
                                                rate_prior_dists       ,
                                                brs              ,
                                                node_path        ,
                                                rate_tuning_par_vec   , rate_pars_acceptfreq_vec, rate_pars_propfreq_vec ,
                                                prior_only
                                                )

        elseif par2upd == :clado_move


        clado_probs_c ,
        edge_probs         ,
        post_clado_probs   ,
        clado_priorlL_c    ,
        prunelL_c          , 
        lL_c               ,
        clado_par_acceptfreq_vec ,
        clado_par_propfreq_vec =clado_par_upd(
                                          rate_pars_c      ,
                                          clado_probs_c     , 
                                          clado_Pmat_c           ,
                                          lL_c                   ,
                                          prunelL_c              ,
                                          clado_priorlL_c        ,
                                          rate_priorlL_c         ,  
                                          edge_probs         ,
                                          post_clado_probs   ,                                               
                                          clado_prior_dists       , 
                                          split_index_vec              ,
                                          sub_index_vec                ,
                                          Q_par_matrix     ,
                                          Q_zeros          ,
                                          Q_index_vec  ,
                                          tip_probs              ,
                                          brs                    ,
                                          node_path              ,
                                          clado_tuning_par_vec   ,
                                          clado_par_acceptfreq_vec     ,
                                          clado_par_propfreq_vec       ,
                                          prior_only
                                          )
        

        elseif par2upd == :rate_zero_switch

            rate_pars_c ,
            rate_pars_switches_c,
            edge_probs  , 
            post_clado_probs,
            rate_priorlL_c   ,
            prunelL_c   , 
            lL_c        ,
            rate_zeroswitch_acceptfreq_vec ,
            rate_zeroswitch_propfreq_vec =rjswitch_rate_par_upd(rate_pars_c      ,
                                                rate_pars_switches_c,
                                                tip_probs        ,
                                                clado_Pmat_c       ,
                                                lL_c             ,
                                                prunelL_c        ,
                                                clado_priorlL_c  ,
                                                rate_priorlL_c   ,  
                                                edge_probs         ,
                                                post_clado_probs   ,                                                
                                                Q_par_matrix     ,
                                                Q_zeros          ,
                                                Q_index_vec  ,
                                                rate_prior_dists       ,
                                                brs              ,
                                                node_path        ,
                                                rate_zeroswitch_acceptfreq_vec, 
                                                rate_zeroswitch_propfreq_vec ,
                                                prior_only
                                                )

        
            
        elseif par2upd == :clado_zero_switch

            clado_probs_c               ,
            clado_probs_switches_c      ,
            edge_probs         ,
            post_clado_probs   ,
            clado_priorlL_c    ,
            prunelL_c          , 
            lL_c               ,
            clado_zeroswitch_acceptfreq_vec ,
            clado_zeroswitch_propfreq_vec =rjswitch_clado_par_upd(rate_pars_c      ,
                                                clado_probs_c      ,
                                                clado_probs_switches_c,
                                                clado_Pmat_c       ,
                                                lL_c                   ,
                                                prunelL_c              ,
                                                clado_priorlL_c        ,
                                                rate_priorlL_c         ,  
                                                edge_probs         ,
                                                post_clado_probs   ,                                               
                                                clado_prior_dists       , 
                                                split_index_vec              ,
                                                sub_index_vec                ,
                                                Q_par_matrix     ,
                                                Q_zeros          ,
                                                Q_index_vec  ,
                                                tip_probs              ,
                                                brs                    ,
                                                node_path              ,
                                                clado_zeroswitch_acceptfreq_vec, 
                                                clado_zeroswitch_propfreq_vec ,
                                                prior_only
                                                )

        
        end


        if mod(iter, anc_state_sampling)==0

            
            anc_states, anc_clado_events=sample_anc_states(edge_probs,  post_clado_probs, root_probs, tip_probs, Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix, clado_Pmat_c, brs, node_path)

            #print(root_probs)
            #print("\n") 

            #print("\n")

            push!(anc_states_cache, join(string.(anc_states),"\t"))
            push!(anc_clado_cache, join(string.(anc_clado_events),"\t"))


            if mod(length(anc_states_cache), write_interval)==0

                anc_states_log_file =  open(anc_states_log_filename,"a")

                [println( anc_states_log_file, anc_states_cache[i]) for i in eachindex(anc_states_cache)]
               
                close( anc_states_log_file )

                anc_states_cache=[]
               
            end


            if mod(length(anc_clado_cache), write_interval)==0

                anc_clado_log_file =  open( anc_clado_log_filename,"a")

                [println( anc_clado_log_file, anc_clado_cache[i]) for i in eachindex(anc_clado_cache)]
               
                close( anc_clado_log_file )

                anc_clado_cache=[]

            end

            
            #[println(log_file, chain_cache[i]) for i in eachindex(chain_cache)]



            #println(anc_clado_log_file, join(string.(anc_clado_events),"\t"))

            #println(anc_state_log_file, join(string.(anc_states),"\t"))
           # print("anc_states logged")
        end


        if mod(iter, iter_trims)==0
            file_line= string(iter) * "\t" *
                       string(round(lL_c, digits=4)) * "\t" * 
                       string(round(prunelL_c,digits=4)) * "\t" * 
                       string(round(rate_priorlL_c,digits=4)) * "\t" * 
                       string(round(clado_priorlL_c, digits=4)) * "\t" *
                       join( string.(round.(rate_pars_c, digits=5)), "\t") * "\t"  * 
                       join( string.(round.(rate_pars_switches_c, digits=5)), "\t") * "\t"  * 

                       join( string.(round.(clado_probs_c, digits=5)), "\t") * "\t"  * 
                       join( string.(round.(clado_probs_switches_c, digits=5)), "\t") * "\t"  * 

                       #string(round(clado_probs_c, digits=3)) * "\t"  * 
                       join( string.(round.(rate_pars_acceptfreq_vec ./rate_pars_propfreq_vec, digits=3)), "\t") * "\t" *
                       join( string.(round.(rate_zeroswitch_acceptfreq_vec ./rate_zeroswitch_propfreq_vec, digits=3)), "\t") * "\t" *

                       join( string.(round.(clado_par_acceptfreq_vec ./clado_par_propfreq_vec, digits=3)), "\t")  * "\t" *
                       join( string.(round.(clado_zeroswitch_acceptfreq_vec ./clado_zeroswitch_propfreq_vec, digits=3)), "\t") 

            
            push!(chain_cache, file_line)


            if mod(length(chain_cache), write_interval)==0

                log_file =  open(log_filename,"a")

                [println(log_file, chain_cache[i]) for i in eachindex(chain_cache)]
               
                close(log_file)
               
                chain_cache=[]
            end
            
            #println(log_file, file_line)
            #print(file_line)
        end

    end 

    end   

    #end
#end
#end




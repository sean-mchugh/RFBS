

function get_job_summary(dir)

    by_doublesingle     = ( sum(occursin.( "_ds",  dir ) ) > 0 )  
    by_gainloss         = ( sum(occursin.( "_gl",  dir ) ) > 0 ) 
    by_rf               = ( sum(occursin.( "_rf",  dir ) ) > 0 )
    by_biome            = ( sum(occursin.( "_b",   dir ) ) > 0 )
    allow_double_gains  = ( sum(occursin.( "2g",   dir ) ) > 0 )
    allow_double_losses = ( sum(occursin.( "2l",   dir ) ) > 0 )
    allow_single_gains  = ( sum(occursin.( "1g",   dir ) ) > 0 )
    allow_single_losses = ( sum(occursin.( "1l",   dir ) ) > 0 )
    DEC                 = ( sum(occursin.( "DEC",   dir ) ) > 0 )  
    ecological          = ( sum(occursin.( "eco",   dir ) ) > 0 )
    allopatric          = ( sum(occursin.( "allo",   dir ) ) > 0 )


    return( by_doublesingle,     
             by_gainloss    ,     
             by_rf          ,     
             by_biome       ,     
             allow_double_gains,  
             allow_double_losses, 
             allow_single_gains , 
             allow_single_losses, 
             DEC                , 
             ecological          ,
             allopatric          )

end

function posterior_tip_pred_wrapper(iters, iter_write, posterior_vec, 
                                    true_tip_states, 
                                    tree_file,
                                    run_name,
                                    model_par_names,
                                    nbiomes,
                                    max_range,
                                    log_filename, append2file, 
                                    job_Bool_vec=[false]
                                    )                        



    #log_header_string= 
    #"samp_state_acc" * "\t" * 
    #"samp_aff_acc" * "\t" *
    #"samp_real_aff_acc" * "\t" * 
    #"samp_fund_aff_acc" * "\t" * 
    #"samp_non_aff_acc" *"\t"*
    #join(model_par_names, "\t")
 
   log_header_string= 
   "r_aff_acc" * "\t" * 
   "total_r_aff_p" * "\t" *
   join([ "r_aff_p_b"*string(i) for i in 1:nbiomes]   , "\t") * "\t" * 
   join(model_par_names, "\t")

    if !isfile(log_filename) 

        log_file =  open(log_filename,"a")
        #write(exampleFileIOStream, file_line)
        println(log_file, log_header_string)
    
        close(log_file)
    
    else 
            if !append2file
                return("ERROR: file already exists and append2file false, cancling post pred")
            end

            print("appending to existing " *log_filename*"\n")
    end
            
    if !(length(job_Bool_vec)==11)

        #performs posterior prediction of tip affinities for a single MCMC 
        #(must provide a RFBS directory file giving the parameters of the MCMC)
        #can skip this and go directly to posterior_tip_pred_sampling if doing this outside of MS pipeline
        by_doublesingle     = ( sum(occursin.( "_ds",  dir ) ) > 0 )  
        by_gainloss         = ( sum(occursin.( "_gl",  dir ) ) > 0 ) 
        by_rf               = ( sum(occursin.( "_rf",  dir ) ) > 0 )
        by_biome            = ( sum(occursin.( "_b",   dir ) ) > 0 )
        allow_double_gains  = ( sum(occursin.( "2g",   dir ) ) > 0 )
        allow_double_losses = ( sum(occursin.( "2l",   dir ) ) > 0 )
        allow_single_gains  = ( sum(occursin.( "1g",   dir ) ) > 0 )
        allow_single_losses = ( sum(occursin.( "1l",   dir ) ) > 0 )
        DEC                 = ( sum(occursin.( "DEC",   dir ) ) > 0 )  
        ecological          = ( sum(occursin.( "eco",   dir ) ) > 0 )
        allopatric          = ( sum(occursin.( "allo",   dir ) ) > 0 )
     
    else

        by_doublesingle     = job_Bool_vec[1]
        by_gainloss         = job_Bool_vec[2]
        by_rf               = job_Bool_vec[3]
        by_biome            = job_Bool_vec[4]
        allow_double_gains  = job_Bool_vec[5]
        allow_double_losses = job_Bool_vec[6]
        allow_single_gains  = job_Bool_vec[7]
        allow_single_losses = job_Bool_vec[8]
        DEC                 = job_Bool_vec[9]
        ecological          = job_Bool_vec[10]
        allopatric          = job_Bool_vec[11]

    end    

    Q_par_matrix, rf_states, move_types, Q_index_vec, move_matrix = make_rf_par_matrix(nbiomes,max_range, 
                                                                                       allow_double_gains , 
                                                                                       allow_double_losses,
                                                                                       allow_single_gains ,
                                                                                       allow_single_losses,
                                                                                       by_biome           ,
                                                                                       by_rf              ,
                                                                                       by_gainloss        ,
                                                                                       by_doublesingle    ,
                                                                                       DEC)

    cladoPmat_unpar=makeclado_Pmat(rf_states, 
                                    ecological,
                                    allopatric)


    String_Clado_Mats ,split_index_vec, sub_index_vec=get_cladoPmat_par_vecs(cladoPmat_unpar)

    #need rate par priors
    posterior_tip_pred_sampling(iters,iter_write,
                                posterior_vec, 
                                true_tip_states, 
                                tree_file,
                                nbiomes,
                                rf_states,
                                model_par_names,
                                Q_par_matrix,  Q_index_vec, 
                                cladoPmat_unpar, split_index_vec, sub_index_vec,
                                log_filename
                                )       

    return("done! written to"*log_filename)                        
end



function posterior_tip_pred_sampling(iters,iter_write,
                                     posterior, 
                                     true_tip_states, 
                                     tree_file,
                                     nbiomes,
                                     rf_states,
                                     model_par_names,
                                     Q_par_matrix,  Q_index_vec, 
                                     cladoPmat_unpar, split_index_vec, sub_index_vec,
                                     log_filename
                                     )                        


             

    
   log_cache=[]

   print("starting simulation"*"\n")

    for iter in 1:iters
         #sample pars from posterior
       #  print(iter)
        posterior_par_sample=posterior[rand(eachindex(posterior))]

        #extract rate pars and clado pars based on par names
        rate_pars_samp=posterior_par_sample[(!).(occursin.( "clado",    model_par_names))]

        clado_pars_samp=posterior_par_sample[(occursin.( "clado",    model_par_names))][1]

        Q_samp=fill_emptyQmat(zeros(size(Q_par_matrix)), rate_pars_samp, Q_index_vec, Q_par_matrix)

        cladoPmat_samp=fill_subsplit_cladoPmat(cladoPmat_unpar, 
                                               clado_pars_samp,
                                               split_index_vec,
                                               sub_index_vec)

        #rtree_file=("R_tree/R_tree_"*string(sample(1:length(readdir("R_tree/")))))
        #rtree_file="R_tree/R_tree_36"

        tip_states_sim, anc_states_sim, clado_events_sim, tree, brs, node_path = sim_rf_tips(rf_states, Q_samp, cladoPmat_samp, tree_file, ntips)

        #samp_state_acc, 
        #samp_aff_acc, 
        #samp_real_aff_acc, 
        #samp_fund_aff_acc,
        #samp_non_aff_acc  =calc_aff_accuracy(true_tip_states, tip_states_sim, rf_states)

        anc_aff_post, 
        anc_aff_post_freq, 
        anc_aff_post_sup,
        anc_aff_acc, 
        aff_correct_mat, 
        true_anc_aff      = calc_anc_state_post_support(rf_states, tip_states_sim', true_tip_states)
    
        real_aff_acc= sum(anc_aff_acc[true_anc_aff .==2])/length(anc_aff_acc[true_anc_aff .==2])
        real_aff_percents= [sum(anc_aff_post_freq[:,i,3]) for i in eachindex(anc_aff_post_freq[1,:,3])]./sum(anc_aff_post_freq[:,:,3])
        total_real_aff_percent= sum(anc_aff_post_freq[:,:,:3])/sum(anc_aff_post_freq[:,:,:])

        log_line_string = string( real_aff_acc)               * "\t" *
                          string(total_real_aff_percent)      * "\t" * 
                          join(real_aff_percents, "\t")       * "\t" *
                          join(posterior_par_sample, "\t")


        #log_line_string= string(samp_state_acc)*"\t"*
        #                 string(samp_aff_acc)*"\t"* 
        #                 string(samp_real_aff_acc)*"\t"*
        #                 string(samp_fund_aff_acc)*"\t"*      
        #                 string(samp_non_aff_acc)*"\t" *
        #                 join(posterior_par_sample, "\t")

        
        
    
        #                 
        push!(log_cache, log_line_string)


        if mod(length(log_cache), iter_write)==0

            log_file =  open(log_filename,"a")

            [println(log_file, log_cache[i]) for i in eachindex(log_cache)]

            close(log_file)

            print(iter)

        end



    end

end 



function calc_aff_accuracy(true_tip_states, tip_states_sim, rf_states)


    true_real=[rf_states[true_tip_states][i].==2 for i in eachindex(true_tip_states)]
    sim_real=[rf_states[tip_states_sim][i].==2 for i in eachindex(tip_states_sim)]
    correct_real_aff=sum([sum( sim_real[i].==true_real[i].==1 ) for i in eachindex(tip_states_sim)])
    real_aff_acc=correct_real_aff/sum(sum.(true_real))


    true_fund=[rf_states[true_tip_states][i].==1 for i in eachindex(true_tip_states)]
    sim_fund=[rf_states[tip_states_sim][i].==1 for i in eachindex(tip_states_sim)]
    correct_fund_aff=sum([sum( sim_fund[i].==true_fund[i].==1 ) for i in eachindex(tip_states_sim)])
    fund_aff_acc=correct_fund_aff/sum(sum.(true_fund))

    if sum(sum.(true_fund))==0
        fund_aff_acc=-9999
    end

    true_non=[rf_states[true_tip_states][i].==0 for i in eachindex(true_tip_states)]
    sim_non=[rf_states[tip_states_sim][i].==0 for i in eachindex(tip_states_sim)]
    correct_non_aff=sum([sum( sim_non[i].==true_non[i].==1 ) for i in eachindex(tip_states_sim)])
    non_aff_acc=correct_non_aff/sum(sum.(true_non))

    if sum(sum.(sum(sum.(true_non))))==0
        non_aff_acc=-9999
    end


    correct_aff=[rf_states[true_tip_states][i].==rf_states[tip_states_sim][i] for i in eachindex(true_tip_states)]
    aff_acc=sum(sum.(correct_aff))/sum(length.(eachindex.(correct_aff)))



    correct_states=true_tip_states.==tip_states_sim
    state_acc=sum(correct_states)/(length(tip_states_sim))

    return(state_acc, aff_acc, real_aff_acc, fund_aff_acc,non_aff_acc )

end


function calc_anc_state_post_support(state_space, anc_post, true_anc_states=Int.(zeros(0)))

    nbiomes=length(state_space[1])
    
    true_anc_aff=Int.(zeros( length(anc_post[1,:]),nbiomes))
    
    anc_aff_acc= zeros( length(anc_post[1,:]),nbiomes)
    
    anc_aff_post_freq=Int.(zeros( length(anc_post[1,:]),nbiomes, nbiomes))
    
    anc_aff_post_sup=zeros( length(anc_post[1,:]),nbiomes, nbiomes)
    
    aff_correct_mat=Int.(zeros(length(anc_post[:,1]), length(anc_post[1,:]),nbiomes))
    
    anc_aff_post=Int.(zeros(length(anc_post[:,1]), length(anc_post[1,:]),nbiomes))
    
    
    for node in eachindex(anc_post[1,:])
         for i in 1:nbiomes 
            for iter in eachindex(anc_post[:,1]) 
    
    
                if (anc_post[iter,node]==0)
                
                else 
                
                    anc_aff_post[iter,node,i]= state_space[anc_post[iter,node]][i] 
                
                    anc_aff_post_freq[node,i,(anc_aff_post[iter,node,i]+1)]+=1
                
                
                    if length(true_anc_states)>0

                        true_anc_aff[node,i]=state_space[true_anc_states[node]][i] 
                
                            if anc_aff_post[iter,node,i]==true_anc_aff[node,i]
                
                                 aff_correct_mat[iter,node,i]=1
                
                            end
        
                    end
                
                end 
            end #end of iter loop


            if (anc_post[1,node]==0)
                
            else 
            
                if length(true_anc_states)>0
                    anc_aff_acc[node,i]=sum(aff_correct_mat[:,node,i])/length(anc_post[:,1])
                end   

                
                anc_aff_post_sup[node,i,:]=anc_aff_post_freq[node,i,:]/length(anc_post[:,1])

            end #end of biome loop


            #print(node)
    
    
        end
    end #end of node loop
    
    
    return(anc_aff_post, anc_aff_post_freq, anc_aff_post_sup,anc_aff_acc, aff_correct_mat, true_anc_aff)
    
end



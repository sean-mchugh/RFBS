########simulated data from prior sampled rates values and estimate with RFBS and DEC

cd(@__DIR__)


#cd("/Volumes/michael.landis/Active/Sean/RFBS/scripts/analysis_scripts")


cd("../..")
pwd()

#ENV["PATH"]="$PATH:/julia-1.9.0/bin"

#uncertain_tips=false
#UP =0.1
#SG =0.1
#SGP=0.1

#print(UP)
#print(SG)
#print(SGP)

#using Pkg; Pkg.add("MKL")

using 
#Tapestree, 
#Tapestree.TRIBE, 
#Tapestree.Utils, 
#Tapestree.ESSE, 
#MKL,
LinearAlgebra,
Pkg, 
Combinatorics , 
StatsBase, 
Distributions, 
Random, 
ExponentialUtilities,
#FastExpm,
# StatsPlots, Plots, 
# Phylo, Optim, 
 DelimitedFiles
 # PProf
#using PhyloNetworks
#using MCMCDiagnosticTools
#using MCMCChains

#Pkg.add("RCall")
#cd("../../../../..")
#ENV["R_HOME"]="/opt/R/4.2.3/lib/R"
#ENV["LD_LIBRARY_PATH"]="`R RHOME`/lib"
#ENV["PATH"]="/opt/R/4.2.3/lib/R"
#ENV["PATH"]="R-4.2.3/bin/R"
#Pkg.build("RCall")

using RCall

include(pwd() *"/src/CladoProbMatrix_fns.jl")
include(pwd() *"/src/make_rf_par_matrix_fns_switch.jl")
include(pwd() *"/src/Tree_Util_fns.jl")
include(pwd() *"/src/rf_sim_fns.jl")
include(pwd() *"/src/rf_rate_clado_mcmc_fns.jl")
include(pwd() *"/src/post_pred_fns.jl")
include(pwd() *"/scripts/scratch_ws/RJ_update_mcmc_ws.jl")

######set run arguements######################

   allow_double_gains =parse(Bool, ARGS[1])
   allow_double_losses=parse(Bool, ARGS[2])
   allow_single_gains =parse(Bool, ARGS[3])
   allow_single_losses=parse(Bool, ARGS[4])
   by_biome           =parse(Bool, ARGS[5])
   by_rf              =parse(Bool, ARGS[6])
   by_gainloss        =parse(Bool, ARGS[7])
   by_doublesingle    =parse(Bool, ARGS[8])
   DEC                =parse(Bool, ARGS[9])                
   UP                 =parse(Float64, ARGS[10])
   SG                 =parse(Float64, ARGS[11])
   SGP                =parse(Float64, ARGS[12])
   uncertain_tips     =parse(Bool,ARGS[13])
   Prior              =parse(Float64,ARGS[14])
   prior_only         =parse(Bool,ARGS[15])
   ntips              =parse(Int64,ARGS[16])
   run_num            =parse(Int64,ARGS[17])           
   
   
   
   #allow_double_gains =true
   #allow_double_losses=false
   #allow_single_gains =true
   #allow_single_losses=true
   #by_biome           =false
   #by_rf              =true
   #by_gainloss        =true
   #by_doublesingle    =true
   #DEC                =false        
   #UP                 =0.1
   #SG                 =0.1
   #SGP                =0.1
   #uncertain_tips     =false
   #Prior=1.0
   #prior_only=false
   #ntips=150
   #run_num=1
   #absorb_states=false

   allow_switches = true

   fun2non=true


    reps=1
    start_truesim=false
    prior_only   

    iters                   =1000000
    iter_trims         =iters/1000
    anc_state_sampling =iters/1000
    write_interval           =1000
    post_pred_iters          =1000
    post_pred_write_interval =1000
    #iters                   =2500
    #iter_trims         =iters/2500
    #anc_state_sampling =iters/2500
    #write_interval           =1
    #post_pred_iters          =1
    #post_pred_write_interval =1


    rate_tuning_par=1.5
    clado_tuning_par=0.4

    N_rate_pars_off              =2
    N_clado_pars_off             =0
    Rtree=true
    #set the value of all pars to 1
    single_par=false

    #iters                   =2500
    #iter_trims         =iters/2500
    #anc_state_sampling =iters/2500
    #write_interval           =1
    #post_pred_iters          =1
    #post_pred_write_interval =1

    proposal_probs = (rate_move = 10 , clado_move = 5 , rate_zero_switch = 5, clado_zero_switch = 0 )


    rate_tuning_par=1.5
    clado_tuning_par=0.4






    clado_types=["sub_split"]

    #DEC=false
    ecological=true
    allopatric=true

    if DEC

        ecological=true
        allopatric=false

    end

    clado_types=["sub_split"]




    #prior label for directory
    prior_label="Exp" * replace(replace(string(Prior),"."=> "p"), "-" => "neg")




######set up priors###############################
    #rate_prior_vec  = [LogNormal(-0.5,0.5)]


    clado_types=["sub_split"]

    #prior_label="LN_neg0p5_0p5"
    prior_label="Exp" * replace(replace(string(Prior),"."=> "p"), "-" => "neg")



######state space and tip uncertainty arguements###############################

    nbiomes_g=[]
    push!(nbiomes_g, 3)

    max_range_g=[]
    push!(max_range_g, nbiomes_g[1])

    #uncertain_tips=true


    uncertain_percent=[]
    down_sample_state_group_percent=[]
    state_group_percent=[]
    uncertain_percent_g    =  [UP]
    state_group_percent_g  =  [SG]
    down_sample_state_group_percent_g  =  [SGP]





    drop_real=true #allow the observed 2's and 0's state to be droppable when nto the true state
    #uncertain_percent_g=Uniform(0,1)



    if uncertain_tips
        for i in 1:reps
            push!(uncertain_percent, rand(uncertain_percent_g))
            push!(state_group_percent, rand(state_group_percent_g))
            push!(down_sample_state_group_percent, rand(down_sample_state_group_percent_g ))

        end
    else 
        for i in 1:reps
            push!(uncertain_percent, 0.0)
            push!(state_group_percent, 0.0)
            push!(down_sample_state_group_percent,  1.0)
        
        end
    end    



######if DEC fix some arguements###############################


    if DEC
        allow_double_gains =true
        allow_double_losses=true
        allow_single_gains =false
        allow_single_losses=false

        #
        by_biome           =false
        by_rf              =false
        by_doublesingle    =false
    end




######Make outfile directory###############################


    dir_name=  string(ntips)*"t_"*string(nbiomes_g[1])*"nB_"*prior_label*"iter_"* string(iters) * "_Ncladooff_" * string(N_clado_pars_off) * "_Nrateoff_" * string(N_rate_pars_off)


    if DEC

        dir_name=dir_name*"_DEC"

    end

    if fun2non
        dir_name=dir_name*"_f2n"
    else
        dir_name=dir_name*"_f2r"

    end

    if start_truesim==true
        dir_name="TS_"*dir_name
    end

    if prior_only==true
        dir_name="PO_"*dir_name
    end

    dir_name=dir_name*"_clado"


    if Rtree==true
        dir_name=dir_name*"_Rtre"
    end
    if allow_double_gains==true
        dir_name=dir_name*"_2g"
    end
    if allow_double_losses==true
        dir_name=dir_name*"_2l"
    end
    if allow_single_gains==true
        dir_name=dir_name*"_1g"
    end

    if allow_single_losses==true
        dir_name=dir_name*"_1l"
    end

    if allow_switches==true
        dir_name=dir_name*"_2sw"
    end


    if single_par==true
        dir_name= dir_name*"_1p"
    end

    if drop_real==true
        dir_name=dir_name*"_dr"
    end


    if uncertain_tips==true
        dir_name=dir_name*"_unce"*string(uncertain_percent_g)*string(state_group_percent_g)*string(down_sample_state_group_percent_g )
    end

    if by_biome==true
        dir_name=dir_name*"_b"
    end
    if by_rf==true
        dir_name= dir_name*"_rf"
    end
    if by_gainloss==true
        dir_name= dir_name*"_gl"
    end


    if by_doublesingle ==true
        dir_name= dir_name*"_ds"
    end


    dir_name ="outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp"*replace(dir_name,r" "  => s"_")


    if !isdir(dir_name)

        if !isdir(dir_name)


            mkpath(dir_name)

        end
    end

######Make model objects########################################################################################


    #make containers for run details

    sim_pars=[]
    estim_pars=[]
    tip_state_sets=[]
    Q_set=[]
    start_sets=[]
    sim_loglik_set=[]
    out_loglik_set=[]
    NaN_run=[]





    print("rep")
    print(" ")
    #print(i)
    print(" ")

    nbiomes=nbiomes_g[1]

    max_range=max_range_g[1]



    allow_real_switches = allow_switches
    allow_realfun_switches = false

    # take details to generater Q vectors to map rates to right Q matrix, the Q matrix generated by this function is ignored (just using the simulation function, messy)
    Q_par_matrix, rf_states, move_types, Q_index_vec, move_matrix = make_rf_par_matrix(nbiomes,max_range,
                                                                                    allow_real_switches     , 
                                                                                    allow_realfun_switches,
                                                                                    allow_double_gains , 
                                                                                    allow_double_losses,
                                                                                    allow_single_gains ,
                                                                                    allow_single_losses,
                                                                                    by_biome           ,
                                                                                    by_rf              ,
                                                                                    by_gainloss        ,
                                                                                    by_doublesingle    ,
                                                                                    DEC )#,absorb_states)



                                                                                    #make range vector to state index dictionaries
    rf2stDict = Dict(copy(rf_states).=>1:length(rf_states))
    st2rfDict = Dict( (1:length(rf_states)).=>copy(rf_states))

    #need rate par priors
    rate_prior_vec   = [Exponential(Prior)]
    rate_prior_dists = make_Prior(rate_prior_vec, move_types)


    #sample rate parameters 
    rate_pars_sim=[round(rand(rate_prior_dists[i]),digits=4) for i in eachindex(move_types)]
    #rate_pars_sim=[0.5883, 0.6475, 5.2363, 1.7844, 1.4718]

   if N_rate_pars_off>0

     off_pars =     sample(eachindex(rate_pars_sim), N_rate_pars_off; replace=false)

     for rate_par in off_pars

         rate_pars_sim[rate_par] = 0.0

     end


   end
   Q_sim=fill_emptyQmat(zeros(size(Q_par_matrix)), rate_pars_sim, Q_index_vec, Q_par_matrix)


    #make cladogenetic probability matrix
    clado_types=["sub_split"]

    #Make base probability matrix
    cladoPmat_unpar= makeclado_Pmat_equal_prob(rf_states,
    ecological,
    allopatric)

    if N_clado_pars_off>0
        off_pars =     sample(1:3, N_clado_pars_off; replace=false)
    else
        off_pars = []
    end
   


    #modify to add additional free parameters, sample "true" cladogenetic paramaters
    cladoPmat_sim, clado_probs_sim, split_index_vec, sub_index_vec, String_Clado_Mats, clado_prior_dists =make_clado_Pmat(rf_states, ecological, allopatric, off_pars)


    file_name=string(run_num)*"__"*join(string.(rate_pars_sim).*"_")*"_"*join(string.(clado_probs_sim[:,2]),"_")



    #get tree from file then simulate data over the tree (including tip data, ancestral states, cladogenetic event types (subset or split) )
   if Rtree==true
       #rtree_file=("R_trees_anc_states_test/R_tree_"*string(ntips))
       if !isdir(dir_name*"/sim_R_trees")
           mkpath(dir_name*"/sim_R_trees")
       end
       rtree_file=("data/sim/R_trees/R_tree_"*string(ntips)*"tips/R_tree_"*string(sample(1:length(readdir("data/sim/R_trees/R_tree_"*string(ntips)*"tips/") ))))
       #rtree_file="R_tree/R_tree_36"
       tip_states, anc_states, clado_events, tree, brs, node_path = sim_rf_tips(rf_states, Q_sim, cladoPmat_sim, rtree_file, ntips)
       #io = open((dir_name*"/"*"Rtrees"*".txt"), "a") 
       #writedlm(io,([file_name*"\t"*join(string.(rtree_file).*"\t")]))                                                                                                                                                                                                                                                                                                                      
       #close(io)  
       io = open((dir_name*"/sim_R_trees/"*file_name*"___"* (split(rtree_file, "/")[2])  * ".txt"), "a") 
       writedlm(io,([file_name*"\t"*join(string.(rtree_file).*"\t")]))                                                                                                                                                                                                                                                                                                                      
       close(io)        
   else 

        tip_states, anc_states, clado_events, tree, brs, node_path = sim_rf_tips(rf_states, Q_sim, cladoPmat_sim, "NA", ntips)

   end


    #convert tip states to state probabilities (each species has a vector of probabilities (n states long)
    tip_probs= get_tip_probs(rf_states,tip_states ,uncertain_tips, uncertain_percent[1], down_sample_state_group_percent[1] ,state_group_percent[1] , drop_real )


######make a bunch of output directories #####################################


    #save anc states from simulation
    if !isdir(dir_name*"/anc_states_sims")
        mkpath(dir_name*"/anc_states_sims")

    end

    if !isdir(dir_name*"/anc_aff_sims")
        mkpath(dir_name*"/anc_aff_sims")

    end


    io = open((dir_name*"/anc_states_sims/"*file_name*"anc_states"*".txt"), "a") 

        writedlm(io,([join(string.(anc_states).*"\t")]))                                                                                                                                                                                                                                                                                                                      
    close(io)        

    #save anc clado events from simulation

    if !isdir(dir_name*"/anc_clados_sims")
        mkpath(dir_name*"/anc_clados_sims")

    end

    io = open((dir_name*"/anc_clados_sims/"*file_name*"anc_clados"*".txt"), "a") 

        writedlm(io,([join(string.(clado_events).*"\t")]))                                                                                                                                                                                                                                                                                                                      
    close(io)        


    #save tip states events from simulation

    if !isdir(dir_name*"/tip_states_sims")
        mkpath(dir_name*"/tip_states_sims")

    end

    io = open((dir_name*"/tip_states_sims/"*file_name*"tip_states"*".txt"), "a") 

        
    writedlm(io,([join(string.(tip_states).*"\t")]))                                                                                                                                                                                                                                                                                                                      
    close(io)        



    if !isfile((dir_name*"/"*"rf_states"*".txt"))


        io = open((dir_name*"/"*"rf_states"*".txt"), "a") 

            [writedlm(io,([i]))     for i in rf_states]                                                                                                                                                                                                                                                                                                                 
        close(io)        

    end




    if !isfile((dir_name*"/"*"run_args"*".txt"))


        io = open((dir_name*"/"*"run_args"*".txt"), "a") 

            writedlm(io, ["allow_double_gains "*string(allow_double_gains )*"\n" *
                     "allow_double_losses"*string(allow_double_losses)*"\n" *
                     "allow_single_gains "*string(allow_single_gains )*"\n" * 
                     "allow_single_losses"*string(allow_single_losses)*"\n" *
                     "by_biome           "*string(by_biome           )*"\n" *           
                     "by_rf              "*string(by_rf              )*"\n" *            
                     "by_gainloss        "*string(by_gainloss        )*"\n" * 
                     "by_doublesingle    "*string(by_doublesingle    )*"\n" *
                     "DEC                "*string(DEC                )*"\n" *
                     "UP                 "*string(UP                 )*"\n" *
                     "SG                 "*string(SG                 )*"\n" *
                     "SGP                "*string(SGP                )*"\n" *
                     "uncertain_tips     "*string(uncertain_tips)]     
                     )                                                                                                                                                                                                                                                                                                                        
        close(io)        

    end





    if !isdir(dir_name*"/logs")
        mkpath(dir_name*"/logs")

    end


######get start pars################################################
    if start_truesim==true

        start_rates =rate_pars_sim
        start_clado_split_prob=clado_split_prob_sim

    else

        start_rates =[rand(rate_prior_dists[i]) for i in eachindex(move_types)]
        #start_clado_split_prob=round.(rand(clado_prior_dists[1]),digits=4) 
        start_clado_probs=sample_sub_split_equal_clado_rates(clado_prior_dists)

    end


    ##

    Q_zeros =zeros(length(rf_states), length(rf_states))
    #file_name=string(rate_pars_sim)




 #   join(string.(clado_events),"\t")

##### run MCMC ######################################


    log_filename=(dir_name*"/logs/"*string(run_num)*"_RFBS_"*file_name*"_log")
    realfun_mcmc(start_rates, start_clado_probs, 
                iters, iter_trims, write_interval,
                anc_state_sampling, 
                Q_par_matrix, Q_index_vec, move_types, rate_prior_dists, 
                cladoPmat_unpar,  clado_prior_dists,
                rf_states,
                tip_probs, 
                tree, 
                log_filename, 
                rate_tuning_par, clado_tuning_par,
                proposal_probs,
                prior_only)



#log_filename=(dir_name*"/"*"RFBS_"*file_name*"_log")
#realfun_mcmc(start_rates, iters, iter_trims,anc_state_sampling, par_matrix, cladoPmat, prior_vec, rf_states, move_types, index_vec, tip_probs, tree, log_filename,tuning_par,prior_only)

# run pruning algorithm to get lielihood for "True" parmaeter values
prunelL_sim, edge_probs,post_clado_probs = rf_prune_algo(Q_sim    ,
tip_probs        ,
cladoPmat_sim       ,
brs              ,
node_path        
         )


     
rate_priorlL_sim = sum([log(pdf(rate_prior_dists[i], rate_pars_sim[i])) for i in eachindex(rate_pars_sim)]) 
clado_priorlL_sim = sum([log(pdf(clado_prior_dists[i],clado_probs_sim[i,2])) for i in eachindex(clado_prior_dists)]) 
#priorlL_c = sum([log(pdf(prior_dists[i], rate_pars_sim[i])) for i in eachindex(rate_pars_sim)]) 

lL_sim= rate_priorlL_sim + clado_priorlL_sim+ prunelL_sim 
#write files
#








burnin=0.5

#
#if !isdir(dir_name*"/post_pred")
#    mkpath(dir_name*"/post_pred")
#end
#

if !isdir(dir_name*"/anc_acc")
    mkpath(dir_name*"/anc_acc")
end

#RFBS_post_pred_log_filename=(dir_name*"/post_pred/"*par_chain_files[1]*"post_pred")

RFBS_post_pred_log_filename=(dir_name*"/post_pred/"*string(run_num)*"_RFBS_"*file_name*"_post_pred")



#RFBS posterior prediction   

RFBS_anc_mat = Int.(readdlm(log_filename*"_anc_states"))
RFBS_log_mat = readdlm(log_filename)
RFBS_par_chain=RFBS_log_mat[Int(round(size(RFBS_log_mat )[1]*burnin)):end,6:(6+length(move_types))]
RFBS_model_par_names=RFBS_log_mat[1,6:(6+length(move_types))]



#anc_aff_post, anc_aff_post_freq,  anc_aff_post_sup, aff_acc_mat, true_anc_aff =calc_anc_state_post_support(rf_states, RFBS_anc_mat, anc_states)



    
    


RFBS_posterior_vec=[Float64.(RFBS_par_chain[i,:]) for i in 1:size(RFBS_par_chain)[1]]


RFBS_job_Bool_vec=[
    by_doublesingle
    by_gainloss 
    by_rf   
    by_biome
    allow_double_gains 
    allow_double_losses
    allow_single_gains 
    allow_single_losses
    DEC
    ecological  
    allopatric  ]



RFBS_anc_states_log_mat = readdlm(dir_name*"/logs/"*string(run_num)*"_RFBS_"*file_name*"_log_anc_states")
RFBS_anc_states_log_mat=Int.(RFBS_anc_states_log_mat)


RFBS_chain_start= Int(round(size(RFBS_anc_states_log_mat)[1]* burnin))
#RFBS_anc_states

RFBS_anc_aff_post, 
RFBS_anc_aff_post_freq, 
RFBS_anc_aff_post_sup,
RFBS_anc_aff_acc, 
RFBS_aff_correct_mat, 
RFBS_true_anc_aff       = calc_anc_state_post_support(rf_states, RFBS_anc_states_log_mat[RFBS_chain_start:end,:], anc_states)
  



#write ancestral biome affinity accuracy to file
log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_RFBS_"*file_name*"_anc_acc","a")

[println(log_file, join(RFBS_anc_aff_acc[i,:],"\t")) for i in eachindex(RFBS_anc_aff_acc[:,1])]

close(log_file)


log_file =  open(dir_name*"/anc_aff_sims/"*string(run_num)*"_RFBS_"*file_name*"_true_anc_aff","a")

[println(log_file, join(RFBS_true_anc_aff[i,:],"\t")) for i in eachindex(RFBS_true_anc_aff[:,1])]

close(log_file)


#write ancestral biome affinity accuracy to file (each group of three rows is post support for (0,1,2) affinity for each biome respectively)
log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_RFBS_"*file_name*"_anc_aff_post_sup","a")

[println(log_file, join(vec(RFBS_anc_aff_post_sup[i,:,:]'),"\t")) for i in eachindex(RFBS_anc_aff_post_sup[:,1,1])]

close(log_file)

###clear object from memory that are very big
RFBS_anc_states_log_mat=nothing                     
RFBS_aff_correct_mat=nothing
RFBS_aff_correct_mat=nothing

                     
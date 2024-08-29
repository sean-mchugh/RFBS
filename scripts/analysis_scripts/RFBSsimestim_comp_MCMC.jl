########simulated data from prior sampled rates values and estimate with RFBS and DEC

#cd(@__DIR__)
LOOCV=false
LOOCV_tip_drop_num=5


 cd(@__DIR__)

cd("..")

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


fund_adj_matrix=[0 1 0;  #biome 1 adjacencies
                 1 0 1;  #biome 2 adjacencies
                 0 1 0]  #biome3  adjacencies

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
Prior        =parse(Float64,ARGS[14])
prior_only   =parse(Bool,ARGS[15])
ntips        =parse(Int64,ARGS[16])

run_num            =parse(Int64,ARGS[17])           
#
#
#
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
fun2non=true


reps=1
start_truesim=false
prior_only   

iters                   =5000000
iter_trims         =iters/5000
anc_state_sampling =iters/5000
write_interval           =1000
post_pred_iters          =5000
post_pred_write_interval =1000
#iters                   =2500
#iter_trims         =iters/2500
#anc_state_sampling =iters/2500
#write_interval           =1
#post_pred_iters          =1
#post_pred_write_interval =1


rate_tuning_par=1.5
clado_tuning_par=0.4





#DEC=false
ecological=true
allopatric=true

if DEC

    ecological=true
    allopatric=false

end


#rate_prior_vec  = [LogNormal(-0.5,0.5)]
rate_prior_vec  = [Exponential(Prior)]


rate_par_proposal_prob=0.8

clado_types=["sub_split"]

#prior_label="LN_neg0p5_0p5"
prior_label="Exp" * replace(replace(string(Prior),"."=> "p"), "-" => "neg")


nbiomes_g=[]
push!(nbiomes_g, 3)

max_range_g=[]
push!(max_range_g, nbiomes_g[1])

Rtree=true
#set the value of all pars to 1
single_par=false
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




#which are non zero Q matrix elements

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





#if absorb_states

#    dir_name="Abs"*string(ntips)*"t_"*string(nbiomes_g[1])*"nB_"*prior_label*"iter_"* string(iters)

#else

    dir_name=string(ntips)*"t_"*string(nbiomes_g[1])*"nB_"*prior_label*"iter_"* string(iters)
#end


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

if LOOCV==true
    dir_name=dir_name*"_LOOCV"
end

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


dir_name ="BSim_RFBS_DEC_comp"*replace(dir_name,r" "  => s"_")


if !isdir(dir_name)
    #getting a file already eixsts error, let this cause a brief random pause before
 for i in 1:rand(10000:100000)

    print("")
    
   
 end   

 if !isdir(dir_name)


 mkdir(dir_name)

 end
end

##########################################################################################################################
#include("/Users/seanmchugh/Projects/realfun_Biome/scripts/RealFunBiomes.jl")

#using .RealFunBiomes




include("CladoProbMatrix_fns.jl")
include("make_rf_par_matrix_fns.jl")
include("Tree_Util_fns.jl")
include("rf_sim_fns.jl")
include("rf_rate_clado_mcmc_fns.jl")
include("post_pred_fns.jl")


#include("CladoProbMatrix_absorb_fns.jl")
#include("make_rf_par_matrix_absorb_fns.jl")
#include("Tree_Util_fns.jl")
#include("rf_sim_absorb_fns.jl")
#include("rf_rate_clado_absorb_mcmc_fns.jl")
#include("post_pred_fns.jl")
##using Profile, ProfileView
#using TimerOutputs
#include("/Users/seanmchugh/Projects/realfun_Biome/scripts/realfun_Biome_fns.jl")




sim_pars=[]
estim_pars=[]
tip_state_sets=[]
Q_set=[]
start_sets=[]
sim_loglik_set=[]
out_loglik_set=[]

NaN_run=[]

par_matrix_g, rf_states_g, move_types_g, index_vec_g, move_matrix_g = make_rf_par_matrix(nbiomes_g[1],nbiomes_g[1], 
                                                                               allow_double_gains , 
                                                                               allow_double_losses,
                                                                               allow_single_gains ,
                                                                               allow_single_losses,
                                                                               by_biome           ,
                                                                               by_rf              ,
                                                                               by_gainloss        ,
                                                                               by_doublesingle    ,
                                                                               DEC)



#tip_state_freq[1][1]




print("rep")
print(" ")
#print(i)
print(" ")

nbiomes=nbiomes_g[1]

max_range=max_range_g[1]




# take details to generater Q vectors to map rates to right Q matrix, the Q matrix generated by this function is ignored (just using the simulation function, messy)
Q_par_matrix, rf_states, move_types, Q_index_vec, move_matrix = make_rf_par_matrix(nbiomes,max_range, 
                                                                                allow_double_gains , 
                                                                                allow_double_losses,
                                                                                allow_single_gains ,
                                                                                allow_single_losses,
                                                                                by_biome           ,
                                                                                by_rf              ,
                                                                                by_gainloss        ,
                                                                                by_doublesingle    ,
                                                                                DEC )#,absorb_states)



rf2stDict = Dict(copy(rf_states).=>1:length(rf_states))
st2rfDict = Dict( (1:length(rf_states)).=>copy(rf_states))

#need rate par priors
rate_prior_dists=make_Prior(rate_prior_vec, move_types)

rate_pars_sim=[round(rand(rate_prior_dists[i]),digits=4) for i in eachindex(move_types)]
#rate_pars_sim=[0.5883, 0.6475, 5.2363, 1.7844, 1.4718]
Q_sim=fill_emptyQmat(zeros(size(Q_par_matrix)), rate_pars_sim, Q_index_vec, Q_par_matrix)


#exp(Q_sim*100)
#Q_sim*sum(brs[:,5])

clado_types=["sub_split"]

cladoPmat_unpar= makeclado_Pmat_equal_prob(rf_states,
ecological,
allopatric)

cladoPmat_sim, clado_probs_sim, split_index_vec, sub_index_vec, String_Clado_Mats, clado_prior_dists =make_clado_Pmat(rf_states, ecological, allopatric )

clado_probs_sim
sum(clado_probs_sim[:,2])

sum(cladoPmat_sim[19,:,:])

[sum(cladoPmat_sim[i,:,:]) for i in 1:length(rf_states)]


file_name=string(run_num)*"__"*join(string.(rate_pars_sim).*"_")*"_"*join(string.(clado_probs_sim[:,2]),"_")

#tip_states, anc_states, tree, brs, node_path = sim_rf_tips(rf_states, Q_sim, cladoPmat, "NA", ntips)

#"0.5883_0.6475_5.2363_1.7844_1.4718_	R_tree/R_tree_36	"

if Rtree==true

    #rtree_file=("R_trees_anc_states_test/R_tree_"*string(ntips))
    if !isdir(dir_name*"/sim_R_trees")
        mkdir(dir_name*"/sim_R_trees")
        
    end

    rtree_file=("R_tree_"*string(ntips)*"tips/R_tree_"*string(sample(1:length(readdir("R_tree_"*string(ntips)*"tips/")))))
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



brs

#tip_states=[7,7,11,5,17,3,19,2,15,15,2,15,15,15,15,15,15,15,15,15,15,2,2,12,18,15,15,15,15,19,15,10,7,7,10,9,9,9,9,4,4,4,9,10,10,7,8,8,8,2,2,2,2,2,14,16,9,9,5,2,4,4,12,18,10,10,10,10,9,9,15,15,15,14,9,9,9,9,9,17,8,8,8,8,8,10,16,7,7,9,1,13,7,7,10,7,16,7,15,9,8,8,9,9,10,5,4,4,4,4,4,4,3,3,3,15,15,8,8,17,14,14,9,9,9,9,9,2,16,15,10,7,10,9,9,4,9,8,10,9,15,8,7,10,17,17,7,7,2,8,9,4,4,4,5,4,4,4,12,4,4,9,9,9,8,8,9,4,2,18,9,10,9,7,9,9,8,9,8,7,1,12,19,10,15,9,9,9,6,6,10,9,9,9,9,9,10,7,4,9,10,8,8,17,8,8,3,3,11,11,11,4,4,4,4,9,17,11,3,7,9,7,10,8,8,8,17,8,8,8,8,8,8,15,17,6,12,4,18,18,12,6,6,3,5,8,8,8,8,17,5,13,1,5,1,5,1,1,1,15,9,9,4,4,8,8,8,14,5,3,14,14,16,14,4,4,4,4,4,2,15,14,3,7,14,17,18,12,4,2,14,16,2,2,16,2,2,1,1,5]

#tip_probs, state_groups,  state_group_tips_vec= get_tip_probs(rf_states,tip_states ,uncertain_tips, uncertain_percent[1], down_sample_state_group_percent[1] ,state_group_percent[1] , drop_real )
#uncertain_tips=true
#uncertain_percent[1]=1.0
#down_sample_state_group_percent[1]=1.0

#state_group_percent[1]=2/3
tip_probs= get_tip_probs(rf_states,tip_states ,uncertain_tips, uncertain_percent[1], down_sample_state_group_percent[1] ,state_group_percent[1] , drop_real )


if LOOCV==true

   LOOCV_tip_ind= sample(eachindex(tip_probs),LOOCV_tip_drop_num, replace=false)
   [tip_probs[LOOCV_tip_ind[i]].=1.0 for i in eachindex(LOOCV_tip_ind)]

end

#[(sum(tip_probs[i])-1)/(length(state_groups[in.(tip_states[i],state_groups)][1])-1) for i in eachindex(tip_probs)]


#test_tip_probs= get_tip_probs(rf_states,tip_states,true, 1.0,1.0, 1.0, drop_real )


#full_uncertain_tip_probs= get_tip_probs(rf_states,tip_states,true, 1.0,1.0, 1.0, drop_real )
#sum(sum.(tip_probs))/sum(sum.(full_uncertain_tip_probs))
#hcat(rf_states,tip_probs[2], ((1:19).==tip_states[2]))

######################


#save anc states from simulation
if !isdir(dir_name*"/anc_states_sims")
    mkdir(dir_name*"/anc_states_sims")
    
end

if !isdir(dir_name*"/anc_aff_sims")
    mkdir(dir_name*"/anc_aff_sims")
    
end


io = open((dir_name*"/anc_states_sims/"*file_name*"anc_states"*".txt"), "a") 

    writedlm(io,([join(string.(anc_states).*"\t")]))                                                                                                                                                                                                                                                                                                                      
close(io)        

#save anc clado events from simulation

if !isdir(dir_name*"/anc_clados_sims")
    mkdir(dir_name*"/anc_clados_sims")
    
end

io = open((dir_name*"/anc_clados_sims/"*file_name*"anc_clados"*".txt"), "a") 

    writedlm(io,([join(string.(clado_events).*"\t")]))                                                                                                                                                                                                                                                                                                                      
close(io)        


#save tip states events from simulation

if !isdir(dir_name*"/tip_states_sims")
    mkdir(dir_name*"/tip_states_sims")
    
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
    mkdir(dir_name*"/logs")
    
end


#get start pars
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




join(string.(clado_events),"\t")

#####run MCMCMC


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
                rate_par_proposal_prob,
                prior_only)



#log_filename=(dir_name*"/"*"RFBS_"*file_name*"_log")
#realfun_mcmc(start_rates, iters, iter_trims,anc_state_sampling, par_matrix, cladoPmat, prior_vec, rf_states, move_types, index_vec, tip_probs, tree, log_filename,tuning_par,prior_only)


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
#    mkdir(dir_name*"/post_pred")
#end
#

if !isdir(dir_name*"/anc_acc")
    mkdir(dir_name*"/anc_acc")
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


  #  anc_aff_post, 
  #  anc_aff_post_freq, 
  #  anc_aff_post_sup,
  #  anc_aff_acc, 
  #  aff_correct_mat, 
  #  true_anc_aff      = calc_anc_state_post_support(rf_states, tip_states', tip_states)
  #  
  #  
  # real_aff_acc= sum(anc_aff_acc[true_anc_aff .==2])/length(anc_aff_acc[true_anc_aff .==2])
  # real_aff_percents= [sum(anc_aff_post_freq[:,i,3]) for i in eachindex(anc_aff_post_freq[1,:,3])]./sum(anc_aff_post_freq[:,:,3])
  # total_real_aff_percent= sum(anc_aff_post_freq[:,:,:3])/sum(anc_aff_post_freq[:,:,:])



#  posterior_tip_pred_wrapper(post_pred_iters, post_pred_write_interval, 
#                           RFBS_posterior_vec, 
#                           tip_states, 
#                           rtree_file,
#                          dir_name,
#                          RFBS_model_par_names,
#                          nbiomes,
#                          max_range,
#                          RFBS_post_pred_log_filename, false,
#                          RFBS_job_Bool_vec
#                          )      
#

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
  



if LOOCV==true

    LOOCV_aff_cc=RFBS_anc_aff_acc[LOOCV_tip_ind,:]

    if !isdir(dir_name*"/LOOCV")
        mkdir(dir_name*"/LOOCV")
    end
    

    log_file =  open(dir_name*"/LOOCV/"*string(run_num)*"_RFBS_"*file_name*"_LOOCV_acc","a")

        [println(log_file, join(LOOCV_aff_cc[i,:],"\t")) for i in eachindex(LOOCV_aff_cc[:,1])]

    close(log_file)


end


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

#copy file used to directory

#cp("scripts/RealFun_script_MCMC.jl", dir_name*"/RealFun_script_MCMC_test"*file_name*".jl")


#end 

#cat([estim_pars[i].-sim_pars[i] for i in eachindex(estim_pars)],dims=2)



#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#######second DEC run ########################################################
#################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################
#
#DEC=true
#ecological=true
#allopatric=true
#
#if DEC
#
#    ecological=true
#    allopatric=false
#
#end
#
#if DEC
#    allow_double_gains =true
#    allow_double_losses=true
#    allow_single_gains =false
#    allow_single_losses=false
#    
#
#    by_biome           =false
#    by_rf              =false
#    by_doublesingle    =false
#end
#
#
#
#DEC_par_matrix, DEC_states, DEC_move_types, DEC_index_vec, DEC_move_matrix = make_rf_par_matrix(nbiomes,max_range, 
#                                                                                allow_double_gains , 
#                                                                                allow_double_losses,
#                                                                                allow_single_gains ,
#                                                                                allow_single_losses,
#                                                                                by_biome           ,
#                                                                                by_rf              ,
#                                                                                by_gainloss        ,
#                                                                                by_doublesingle    ,
#                                                                                DEC)
##
#
#
#
#
#
#rf_states
#
#rf_states
#
#
#
#
#rf_state_vec=tip_states
#
#DEC_tip_states=convert_RF_2_DEC_states(tip_states,rf_states, DEC_states)
#
#DEC_anc_states=convert_RF_2_DEC_states(anc_states, rf_states, DEC_states)
#
#
#DEC_tip_probs= get_tip_probs(DEC_states,DEC_tip_states,uncertain_tips, uncertain_percent[1],state_group_percent[1], state_group_percent[1], drop_real )
#
#
#DEC_start_rates =[rand(rate_prior_dists[i]) for i in eachindex(DEC_move_types)]
#DEC_start_clado_split_prob=round(rand(clado_prior_dists[1]),digits=4) 
#
#
##DEC_start_rates =[rand(rate_prior_dists[i]) for i in eachindex(DEC_move_types)]
#
##DEC_prior_dists=make_Prior(prior_vec, DEC_move_types)
#
#DEC_cladoPmat_unpar=makeclado_Pmat(DEC_states, ecological, allopatric)
#
#String_Clado_Mats ,DEC_split_index_vec, DEC_sub_index_vec=get_cladoPmat_par_vecs(DEC_cladoPmat_unpar)
#
#cladoPmat_sim=fill_subsplit_cladoPmat(DEC_cladoPmat_unpar, 
#                                      DEC_start_clado_split_prob,
#                                      DEC_split_index_vec,
#                                      DEC_sub_index_vec)
#
#
##log_filename=(dir_name*"/"*"DEC_"*file_name*"_log")
##realfun_mcmc(DEC_start_rates, iters, iter_trims,
##            anc_state_sampling, DEC_par_matrix, 
##            DEC_cladoPmat, prior_vec, DEC_states, DEC_move_types, 
##           index_vec, DEC_tip_probs, tree, 
##            log_filename,tuning_par,prior_only)
##
#DEC_log_filename=(dir_name*"/logs/"*string(run_num)*"_DEC_"*file_name*"_log")
#realfun_mcmc(DEC_start_rates, start_clado_split_prob, 
#            iters, iter_trims, write_interval,
#            anc_state_sampling, 
#            DEC_par_matrix, DEC_index_vec,  DEC_move_types, rate_prior_vec, 
#            DEC_cladoPmat_unpar, clado_types, clado_prior_vec,
#            DEC_states,
#            DEC_tip_probs, 
#            tree, 
#            DEC_log_filename, 
#            rate_tuning_par, clado_tuning_par,
#            rate_par_proposal_prob,
#            prior_only)
#
#
#
#
#
#
#
#if !isdir(dir_name*"/post_pred")
#    mkdir(dir_name*"/post_pred")
#end
#
#
#DEC_post_pred_log_filename=(dir_name*"/post_pred/"*string(run_num)*"_DEC_"*file_name*"_post_pred")
#
#
##DEC posterior prediction  
#DEC_log_mat = readdlm(DEC_log_filename)
#DEC_par_chain=DEC_log_mat[Int(round(size(DEC_log_mat )[1]*burnin)):end,6:(6+length(DEC_move_types))]
#DEC_model_par_names=DEC_log_mat[1,6:(6+length(DEC_move_types))]
#
#
#DEC_posterior_vec=[Float64.(DEC_par_chain[i,:]) for i in 1:size(DEC_par_chain)[1]]
#  
#
#DEC_job_Bool_vec=[
#    by_doublesingle
#    by_gainloss 
#    by_rf   
#    by_biome
#    allow_double_gains 
#    allow_double_losses
#    allow_single_gains 
#    allow_single_losses
#    DEC
#    ecological  
#    allopatric  ]
#    
#
#posterior_tip_pred_wrapper(post_pred_iters,post_pred_write_interval, DEC_posterior_vec, 
#                            DEC_tip_states, 
#                            rtree_file,
#                           dir_name,
#                           DEC_model_par_names,
#                           nbiomes,
#                           max_range,
#                           DEC_post_pred_log_filename, false,
#                           DEC_job_Bool_vec
#                           )      
#
#
#DEC_anc_states_log_mat = readdlm(dir_name *"/logs/"*string(run_num)* "_DEC_"*file_name*"_log_anc_states")
#
#DEC_anc_states_log_mat=Int.(DEC_anc_states_log_mat)
#DEC_chain_start= Int(round(size(DEC_anc_states_log_mat)[1]* burnin))
#
#DEC_anc_states_true=convert_RF_2_DEC_states(anc_states, rf_states, DEC_states)#[(ntips+1):end]
#                       
#
#
#DEC_anc_aff_post, 
#DEC_anc_aff_post_freq, 
#DEC_anc_aff_post_sup,
#DEC_anc_aff_acc, 
#DEC_aff_correct_mat, 
#DEC_true_anc_aff = calc_anc_state_post_support(DEC_states, DEC_anc_states_log_mat[DEC_chain_start:end,:], DEC_anc_states)
#
##write ancestral biome affinity accuracy to file
#log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_DEC_"*file_name*"_anc_acc","a")
#
#[println(log_file, join(DEC_anc_aff_acc[i,:],"\t")) for i in eachindex(DEC_anc_aff_acc[:,1])]
#
#close(log_file)
#
##write ancestral biome affinity accuracy to file
#log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_DEC_"*file_name*"_anc_aff_post_sup","a")
#
#[println(log_file, join(vec(DEC_anc_aff_post_sup[i,:,:]'),"\t")) for i in eachindex(DEC_anc_aff_post_sup[:,1,1])]
#
#close(log_file)
#
#
#
#
#if LOOCV==true
#
#    LOOCV_aff_cc=DEC_anc_aff_acc[LOOCV_tip_ind,:]
#
#    if !isdir(dir_name*"/LOOCV")
#        mkdir(dir_name*"/LOOCV")
#    end
#    
#
#    log_file =  open(dir_name*"/LOOCV/"*"DEC_"*file_name*"_LOOCV_acc","a")
#
#        [println(log_file, join(LOOCV_aff_cc[i,:],"\t")) for i in eachindex(LOOCV_aff_cc[:,1])]
#
#    close(log_file)
#
#
#end
#
#
##DEC_post_pred_log_filename=(dir_name*"/post_pred/"*"DEC_"*file_name*"_post_pred")
###DEC posterior prediction            
##DEC_log_mat = readdlm(DEC_log_filename)
##DEC_par_chain=DEC_log_mat[Int(round(size(DEC_log_mat )[1]*burnin)):end,6:(6+length(DEC_move_types))]
##DEC_model_par_names=DEC_log_mat[1,6:(6+length(DEC_move_types))]
##
##DEC_posterior_vec=[Float64.(DEC_par_chain[i,:]) for i in 1:size(DEC_par_chain)[1]]
                        
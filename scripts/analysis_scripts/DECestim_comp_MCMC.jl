
fund_adj_matrix=[0 1 0;  #biome 1 adjacencies
1 0 1;  #biome 2 adjacencies
0 1 0]  #biome3  adjacencies



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
fun2non=true


#absorb_states=false





#cd(@__DIR__)
LOOCV=false
LOOCV_tip_drop_num=5


cd(@__DIR__)

cd("..")

#uncertain_tips=false
#UP =0.1
#SG =0.1
#SGP=0.1
 # PProf
#using PhyloNetworks
#using MCMCDiagnosticTools
#using MCMCChains

reps=1
start_truesim=false


iters                  =1000000
iter_trims         =iters/5000
anc_state_sampling =iters/5000
write_interval           =1000
post_pred_iters          =5000
post_pred_write_interval =1000

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


#clado_prior_vec = [Beta(1,1)]

rate_par_proposal_prob=0.8

clado_types=["sub_split"]

#prior_label="LN_neg0p5_0p5"
prior_label="Exp" * replace(string(Prior),"."=> "p")



#prior_label="LN_neg0p5_0p5"

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
#
#    dir_name="Abs"*string(ntips)*"t_"*string(nbiomes_g[1])*"nB_"*prior_label*"iter_"* string(iters)
#
#else
#
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



    ecological=true
    allopatric=false


    allow_double_gains =true
    allow_double_losses=true
    allow_single_gains =false
    allow_single_losses=false
    

    by_biome           =false
    by_rf              =false
    by_doublesingle    =false


#####################################################################################################

#include("CladoProbMatrix_absorb_fns.jl")
#include("make_rf_par_matrix_absorb_fns.jl")
#include("Tree_Util_fns.jl")
#include("rf_sim_absorb_fns.jl")
#include("rf_rate_clado_absorb_mcmc_fns.jl")
#include("post_pred_fns.jl")
include("CladoProbMatrix_fns.jl")
include("make_rf_par_matrix_fns.jl")
include("Tree_Util_fns.jl")
include("rf_sim_fns.jl")
include("rf_rate_clado_mcmc_fns.jl")
include("post_pred_fns.jl")



####get phylogeny##

simtree_files=readdir(dir_name*"/sim_R_trees/")



rtree_label_file=simtree_files[startswith.( simtree_files, string(run_num)*"_")][1]


rtree_file=("R_tree_"*string(ntips)*"tips/R_tree_"*split(split(simtree_files[startswith.( simtree_files, string(run_num))][1], "_")[end], ".")[1] )

tree, bts = read_tree(rtree_file)  



file_name=split(rtree_label_file, "___")[1]

nbiomes=3
max_range=3


DEC_par_matrix, DEC_states, DEC_move_types, DEC_index_vec, DEC_move_matrix = make_rf_par_matrix(nbiomes,max_range, 
                                                                                allow_double_gains , 
                                                                                allow_double_losses,
                                                                                allow_single_gains ,
                                                                                allow_single_losses,
                                                                                by_biome           ,
                                                                                by_rf              ,
                                                                                by_gainloss        ,
                                                                                by_doublesingle    ,
                                                                                true) #,absorb_states)
#



rate_prior_dists=make_Prior(rate_prior_vec, DEC_move_types)
#clado_prior_dists=make_Prior(clado_prior_vec, clado_types)




#RFBS_states=make_rf_states(nbiomes,max_range, false, absorb_states)

RFBS_states=make_rf_states(nbiomes,max_range, false)

tip_states_files=readdir(dir_name*"/tip_states_sims")
tip_states=parse.(Int64,split(readdlm(dir_name*"/tip_states_sims/" * tip_states_files[startswith.(tip_states_files, string(run_num))][1])[1], "\t")[1:(ntips)])

anc_states_files=readdir(dir_name*"/anc_states_sims")

RFBS_anc_states=parse.(Int64,split( readdlm( dir_name*"/anc_states_sims/" * anc_states_files[startswith.(anc_states_files, string(run_num))][1] )[1], "\t" )[1:(ntips*4)-2] )


rf_state_vec=tip_states


DEC_tip_states=convert_RF_2_DEC_states(tip_states,RFBS_states, DEC_states, fun2non)
#DEC_anc_states=convert_RF_2_DEC_states(RFBS_anc_states, RFBS_states, DEC_states)


#DEC_anc_states=convert_RF_2_DEC_states(anc_states, rf_states, DEC_states)


DEC_tip_probs= get_tip_probs(DEC_states,DEC_tip_states,uncertain_tips, uncertain_percent[1],state_group_percent[1], state_group_percent[1], drop_real )


DEC_start_rates =[rand(rate_prior_dists[i]) for i in eachindex(DEC_move_types)]
#DEC_start_clado_split_prob=round(rand(clado_prior_dists[1]),digits=4) 


#DEC_start_rates =[rand(rate_prior_dists[i]) for i in eachindex(DEC_move_types)]

#DEC_prior_dists=make_Prior(prior_vec, DEC_move_types)

DEC_cladoPmat_unpar=makeclado_Pmat_equal_prob(DEC_states, ecological, allopatric)

#String_Clado_Mats ,DEC_split_index_vec, DEC_sub_index_vec=get_cladoPmat_par_vecs(DEC_cladoPmat_unpar)

#cladoPmat_sim=fill_subsplit_cladoPmat(DEC_cladoPmat_unpar, 
#                                      DEC_start_clado_split_prob,
#                                      DEC_split_index_vec,
#                                      DEC_sub_index_vec)
#
#                                      clado_types=["sub_split"]

#                                      DEC_cladoPmat_unpar= makeclado_Pmat_equal_prob(rf_states,
#ecological,
#allopatric)

cladoPmat_start, DEC_clado_probs_start, DEC_split_index_vec, DEC_sub_index_vec, String_Clado_Mats, clado_prior_dists =make_clado_Pmat(DEC_states, ecological, allopatric )


#log_filename=(dir_name*"/"*"DEC_"*file_name*"_log")
#realfun_mcmc(DEC_start_rates, iters, iter_trims,
#            anc_state_sampling, DEC_par_matrix, 
#            DEC_cladoPmat, prior_vec, DEC_states, DEC_move_types, 
#           index_vec, DEC_tip_probs, tree, 
#            log_filename,tuning_par,prior_only)
#
DEC_log_filename=(dir_name*"/logs/"*string(run_num)*"_DEC_"*file_name*"_log")
realfun_mcmc(DEC_start_rates,DEC_clado_probs_start, 
            iters, iter_trims, write_interval,
            anc_state_sampling, 
            DEC_par_matrix, DEC_index_vec,  DEC_move_types, rate_prior_dists, 
            DEC_cladoPmat_unpar, clado_prior_dists,
            DEC_states,
            DEC_tip_probs, 
            tree, 
            DEC_log_filename, 
            rate_tuning_par, clado_tuning_par,
            rate_par_proposal_prob,
            prior_only)





#
#
#if !isdir(dir_name*"/post_pred")
#    mkdir(dir_name*"/post_pred")
#end
#

burnin=0.5

DEC_post_pred_log_filename=(dir_name*"/post_pred/"*string(run_num)*"_DEC_"*file_name*"_post_pred")


#DEC posterior prediction  
DEC_log_mat = readdlm(DEC_log_filename)
DEC_par_chain=DEC_log_mat[Int(round(size(DEC_log_mat )[1]*burnin)):end,6:(6+length(DEC_move_types))]
DEC_model_par_names=DEC_log_mat[1,6:(6+length(DEC_move_types))]


DEC_posterior_vec=[Float64.(DEC_par_chain[i,:]) for i in 1:size(DEC_par_chain)[1]]
  

DEC_job_Bool_vec=[
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


DEC_anc_states_log_mat = readdlm(dir_name *"/logs/"*string(run_num)* "_DEC_"*file_name*"_log_anc_states")

DEC_anc_states_log_mat=Int.(DEC_anc_states_log_mat)
DEC_chain_start= Int(round(size(DEC_anc_states_log_mat)[1]* burnin))

DEC_anc_states_true=convert_RF_2_DEC_states(RFBS_anc_states, RFBS_states, DEC_states)#[(ntips+1):end]
                       


DEC_anc_aff_post, 
DEC_anc_aff_post_freq, 
DEC_anc_aff_post_sup,
DEC_anc_aff_acc, 
DEC_aff_correct_mat, 
DEC_true_anc_aff = calc_anc_state_post_support(DEC_states, DEC_anc_states_log_mat[DEC_chain_start:end,:], DEC_anc_states_true)

#write ancestral biome affinity accuracy to file
log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_DEC_"*file_name*"_anc_acc","a")

[println(log_file, join(DEC_anc_aff_acc[i,:],"\t")) for i in eachindex(DEC_anc_aff_acc[:,1])]

close(log_file)

#write ancestral biome affinity accuracy to file
log_file =  open(dir_name*"/anc_acc/"*string(run_num)*"_DEC_"*file_name*"_anc_aff_post_sup","a")

[println(log_file, join(vec(DEC_anc_aff_post_sup[i,:,:]'),"\t")) for i in eachindex(DEC_anc_aff_post_sup[:,1,1])]

close(log_file)




if LOOCV==true

    LOOCV_aff_cc=DEC_anc_aff_acc[LOOCV_tip_ind,:]

    if !isdir(dir_name*"/LOOCV")
        mkdir(dir_name*"/LOOCV")
    end
    

    log_file =  open(dir_name*"/LOOCV/"*"DEC_"*file_name*"_LOOCV_acc","a")

        [println(log_file, join(LOOCV_aff_cc[i,:],"\t")) for i in eachindex(LOOCV_aff_cc[:,1])]

    close(log_file)


end


#DEC_post_pred_log_filename=(dir_name*"/post_pred/"*"DEC_"*file_name*"_post_pred")
##DEC posterior prediction            
#DEC_log_mat = readdlm(DEC_log_filename)
#DEC_par_chain=DEC_log_mat[Int(round(size(DEC_log_mat )[1]*burnin)):end,6:(6+length(DEC_move_types))]
#DEC_model_par_names=DEC_log_mat[1,6:(6+length(DEC_move_types))]
#
#DEC_posterior_vec=[Float64.(DEC_par_chain[i,:]) for i in 1:size(DEC_par_chain)[1]]
                        
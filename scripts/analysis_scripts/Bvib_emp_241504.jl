#cd(@__DIR__)

#cd("..")





#which are non zero Q matrix elements
allow_double_gains =parse(Bool, ARGS[1])
allow_double_losses=parse(Bool, ARGS[2])
allow_single_gains =parse(Bool, ARGS[3])
allow_single_losses=parse(Bool, ARGS[4])
by_biome           =parse(Bool, ARGS[5])
by_rf              =parse(Bool, ARGS[6])
by_gainloss        =parse(Bool, ARGS[7])
by_doublesingle    =parse(Bool, ARGS[8])
DEC=parse(Bool, ARGS[9])
include_fund=parse(Bool, ARGS[10])

use_forbidden_fund_vec=parse(Bool, ARGS[11])
use_adj_matrix=parse(Bool, ARGS[12])
biome_exc_treatment_ind = parse(Int, ARGS[13])
biome_inc_treatment_ind =  parse(Int, ARGS[14])
run_num=parse(Int, ARGS[15])


#allow_double_gains =true
#allow_double_losses=false
#allow_single_gains =true
#allow_single_losses=true
#by_biome           =false
#by_rf              =false
#by_gainloss        =true
#by_doublesingle    =true
#DEC                =false     
#UP                 =0.1
#SG                 =0.1
#SGP                =0.1
#uncertain_tips     =true
#use_adj_matrix=true
#use_forbidden_fund_vec=true
#include_fund=true
#biome_exc_treatment_ind=1
#biome_exc_treatment_ind=2
#run_num=1

using 
#Tapestree, 
#Tapestree.TRIBE, 
#Tapestree.Utils, 
#Tapestree.ESSE, 
Pkg, 
Combinatorics , 
StatsBase, 
Distributions, 
Random, RCall,
ExponentialUtilities,
# StatsPlots, Plots, 
# Phylo, Optim, 
 DelimitedFiles
 # PProf
#using MCMCDiagnosticTools
#using MCMCChains


reps=1
#start_truesim=true
Fossil=true

biome_inc_treatment = ["cons" "bold"][biome_inc_treatment_ind]

if include_fund


   # viburnum_dat=readdlm("viburnum_data_files/incf_viburnum_sorted_rf_states_3b.txt", Int64)
   viburnum_dat=readdlm("viburnum_data_files/"*biome_inc_treatment*"_incf_viburnum_sorted_rf_states_3b.txt", Int64)

    if Fossil
        tree_file="viburnum_data_files/out.1.t163.f5.mask_fossil_states.mcc.tre" #pollen fossil taxa
    else

        tree_file="viburnum_data_files/viburnum_sorted.tre" #extant only
    end
    tip_rf_ranges=[viburnum_dat[i,:] for i in 1:length(viburnum_dat[:,1])]


else

    viburnum_dat=readdlm("viburnum_data_files/viburnum_sorted_rf_states_3b.txt", Int64)
    if Fossil
        tree_file="viburnum_data_files/out.1.t163.f5.mask_fossil_states.mcc.tre" #pollen fossil taxa
    else

        tree_file="viburnum_data_files/viburnum_sorted.tre" #extant only
    end
    tip_rf_ranges=[viburnum_dat[i,:] for i in 1:length(viburnum_dat[:,1])]
end



biome_exc_treatment_vec = [
#"3.biomes.germination.only"                                      
#"3.biomes.germination.bold"                                      
#"3.biomes.germination.conservative"                              
#"3.biomes.leafing.bold"                                          
#"3.biomes.leafing.conservative"                                  
#"3.biomes.USDA"                                                  
#"3.biomes.germination.only.germination.bold"                     
#"3.biomes.germination.only.leafing.bold"                         
#"3.biomes.germination.only.leafing.conservative"                 
#"3.biomes.germination.only.USDA"                                 
#"3.biomes.germination.bold.germination.conservative"             
#"3.biomes.germination.bold.leafing.bold"                         
#"3.biomes.germination.bold.leafing.conservative"                 
#"3.biomes.germination.bold.USDA"                                 
#"3.biomes.germination.conservative.leafing.bold"                 
#"3.biomes.germination.conservative.leafing.conservative"         
#"3.biomes.germination.conservative.USDA"                         
#"3.biomes.leafing.bold.leafing.conservative"                     
#"3.biomes.leafing.bold.USDA"                                     
#"3.biomes.leafing.conservative.USDA"                             
"3.biomes.germination.only.germination.bold.leafing.bold"        
#"3.biomes.germination.only.germination.bold.leafing.conservative"
#"3.biomes.germination.only.germination.bold.USDA"                
#"3.biomes.germination.only.leafing.bold.leafing.conservative"    
"3.biomes.germination.only.leafing.conservative.USDA"            
#"3.biomes.germination.bold.germination.conservative.USDA"        
#"3.biomes.germination.conservative.leafing.conservative.USDA"    
#"3.biomes.leafing.bold.leafing.conservative.USDA"           
]


biome_exc_treatment = biome_exc_treatment_vec[biome_exc_treatment_ind]


uncertain_tips=true
drop_real=false

parameterized_clado=true

uncertain_percent=[]
state_group_percent=[]

uncertain_percent_g    =  [1.0]
state_group_percent_g  =  [1.0]


nbiomes_g=[]
push!(nbiomes_g, length(tip_rf_ranges[1]))
max_range_g=[]
push!(max_range_g, nbiomes_g[1])


if use_adj_matrix

    if nbiomes_g[1]==3

        fund_adj_matrix=[0 1 0;  #biome 1 adjacencies
                         1 0 1 ;  #biome 2 adjacencies
                         0 1 0 ]  #biome3  adjacencies

    end

else
    fund_adj_matrix=NaN
end


if use_forbidden_fund_vec

    
    if nbiomes_g[1]==3
    
         forbidden_aff = Int.(readdlm("viburnum_data_files/"*biome_exc_treatment*".csv", ',', header=true)[1])
    
 
    end

else
    forbidden_aff=NaN
end


iters                 =2000000
iter_trims        =iters/10000
anc_state_sampling=iters/10000
write_interval           =1000
post_pred_iters         =10000
post_pred_write_interval =1000


rate_tuning_par=1.5
clado_tuning_par=0.4


#ntips=10

ecological=true
allopatric=true

if DEC

    ecological=true
    allopatric=false

end

Prior=0.5
rate_prior_vec  = [Exponential(Prior)]
#clado_prior_vec = [Beta(1,1)]

rate_par_proposal_prob=0.8

clado_types=["sub_split"]


prior_label="Exp" * replace(replace(string(Prior),"."=> "p"), "-" => "neg")

prior_only=false
#Rtree=true
#set the value of all pars to 1
#single_par=false
#uncertain_percent_g=Uniform(0,1)





if uncertain_tips
    for i in 1:reps
        push!(uncertain_percent  , rand(uncertain_percent_g))
        push!(state_group_percent, rand(state_group_percent_g))

    end
else 
    for i in 1:reps
        push!(uncertain_percent, 0.0)
        push!(state_group_percent, 1.0)
    end
end    





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



#if Rtree==true
#
#    ntips=300
#
#end


dir_name=string(nbiomes_g[1])*"nB_"*prior_label*"_"* string(iters)


if prior_only==true
    dir_name=dir_name*"_PO"
end

if Fossil==true
    dir_name=dir_name*"_Foss"
end


#if start_truesim==true
#    dir_name="TS_"*dir_name
#end

if use_adj_matrix==true
    dir_name=dir_name*"_admat"
end

if include_fund==true
    #dir_name=dir_name*"_incf"
    dir_name=dir_name*"_incf"
    dir_name=dir_name * "_" * biome_inc_treatment
    


end


if use_forbidden_fund_vec==true
    dir_name=dir_name*"_excf"
    dir_name=dir_name * "_" * biome_exc_treatment

end



if DEC==true
    dir_name=dir_name*"_DEC"
end
if ecological==true
    dir_name=dir_name*"_eco"
end
if allopatric==true
    dir_name=dir_name*"_allo"
end

dir_name=dir_name*"_clado"


#if Rtree==true
#    dir_name=dir_name*"_Rtre"
#end
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
#if single_par==true
#    dir_name= dir_name*"_1p"
#end
#if uncertain_tips==true
#    dir_name=dir_name*"_unce"*string(uncertain_percent_g)*string(state_group_percent_g)
#end
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


if Fossil==true
    dir_name=dir_name*"_Foss"
end

dir_name ="Bvib_"*replace(dir_name,r" "  => s"_")


if !isdir(dir_name)
    mkdir(dir_name)
    
end


#include("/Users/seanmchugh/Projects/realfun_Biome/scripts/RealFunBiomes.jl")

#using .RealFunBiomes



include("CladoProbMatrix_fns.jl")
include("make_rf_par_matrix_fns.jl")
include("Tree_Util_fns.jl")
include("rf_sim_fns.jl")
include("rf_rate_clado_mcmc_fns.jl")
include("post_pred_fns.jl")
#using Profile, ProfileView
#using TimerOutputs
#include("/Users/seanmchugh/Projects/realfun_Biome/scripts/realfun_Biome_fns.jl")



##


#sim_pars=[]

estim_pars=[]
#tip_state_sets=[]
Q_set=[]
start_sets=[]
#sim_loglik_set=[]
out_loglik_set=[]

NaN_run=[]



#
print("rep")
print(" ")
#print(i)
print(" ")

nbiomes=nbiomes_g[1]

max_range=max_range_g[1]

 # take details to generater Q fectors to map rates to right Q matrix, the Q matrix generated by this function is ignored (just using the simulation function, messy)
par_matrix, rf_states, move_types, index_vec, move_matrix = make_rf_par_matrix(nbiomes,max_range, 
                                                                                allow_double_gains , 
                                                                                allow_double_losses,
                                                                                allow_single_gains ,
                                                                                allow_single_losses,
                                                                                by_biome           ,
                                                                                by_rf              ,
                                                                                by_gainloss        ,
                                                                                by_doublesingle    ,
                                                                                DEC)




#need rate par priors

print(tip_rf_ranges)

rf2stDict = Dict(copy(rf_states).=>1:length(rf_states))
st2rfDict = Dict( (1:length(rf_states)).=>copy(rf_states))

tip_states=[rf2stDict[tip_rf_ranges[i]] for i in eachindex(tip_rf_ranges)]

cladoPmat_unpar= makeclado_Pmat_equal_prob(rf_states,
ecological,
allopatric)


file_name="Bvib__"*string(run_num)*"_"


cladoPmat_sim, clado_probs_start, split_index_vec, sub_index_vec, String_Clado_Mats, clado_prior_dists =make_clado_Pmat(rf_states, ecological, allopatric )


rate_prior_dists=make_Prior(rate_prior_vec, move_types)


#clado_types=["sub_split"]
#clado_prior_dists=make_Prior(clado_prior_vec, clado_types)



tree, bts = read_tree(tree_file) 
ntips=length(tree.tlab)
ed=tree.ed
edges = cat(ed, [2*ntips ntips + 1], dims = 1)
triads=maketriads(edges)
node_path  = get_trav_path_4prun(triads )
br = branching_times(tree)
brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
 



#if Rtree==true
#
#    #rtree_file=("R_trees_anc_states_test/R_tree_"*string(ntips))
#
#    rtree_file=("R_tree/R_tree_"*string(sample(1:length(readdir("R_tree/")))))
#    #rtree_file="R_tree/R_tree_36"
#
#    tip_states, anc_states, tree, brs, node_path = sim_rf_tips(rf_states, Q_sim, cladoPmat, rtree_file, ntips)
#
#    io = open((dir_name*"/"*"Rtrees"*".txt"), "a") 
#    writedlm(io,([file_name*"\t"*join(string.(rtree_file).*"\t")]))                                                                                                                                                                                                                                                                                                                      
#    close(io)        
#    
#
#else 
#
#    tip_states, anc_states, tree, brs, node_path = sim_rf_tips(rf_states, Q_sim, cladoPmat, "NA", ntips)
#
#end

#tip_states=[7,7,11,5,17,3,19,2,15,15,2,15,15,15,15,15,15,15,15,15,15,2,2,12,18,15,15,15,15,19,15,10,7,7,10,9,9,9,9,4,4,4,9,10,10,7,8,8,8,2,2,2,2,2,14,16,9,9,5,2,4,4,12,18,10,10,10,10,9,9,15,15,15,14,9,9,9,9,9,17,8,8,8,8,8,10,16,7,7,9,1,13,7,7,10,7,16,7,15,9,8,8,9,9,10,5,4,4,4,4,4,4,3,3,3,15,15,8,8,17,14,14,9,9,9,9,9,2,16,15,10,7,10,9,9,4,9,8,10,9,15,8,7,10,17,17,7,7,2,8,9,4,4,4,5,4,4,4,12,4,4,9,9,9,8,8,9,4,2,18,9,10,9,7,9,9,8,9,8,7,1,12,19,10,15,9,9,9,6,6,10,9,9,9,9,9,10,7,4,9,10,8,8,17,8,8,3,3,11,11,11,4,4,4,4,9,17,11,3,7,9,7,10,8,8,8,17,8,8,8,8,8,8,15,17,6,12,4,18,18,12,6,6,3,5,8,8,8,8,17,5,13,1,5,1,5,1,1,1,15,9,9,4,4,8,8,8,14,5,3,14,14,16,14,4,4,4,4,4,2,15,14,3,7,14,17,18,12,4,2,14,16,2,2,16,2,2,1,1,5]

#tip_probs= get_tip_probs(rf_states, tip_states, uncertain_tips, uncertain_percent[1],state_group_percent[1], drop_real )

tip_probs = get_emp_tip_probs(rf_states, tip_states, fund_adj_matrix, forbidden_aff)


#rf_states[tip_states[158]]
#forbidden_aff[158,:]

#rf_states[tip_probs[158].==1.0]

#rf_states[tip_probs[tree.tlab.=="V_clemensiae"][1].==1.0]
#rf_states[edge_probs[161][1].==1.0]


tree.tlab
rf_states
full_uncertain_tip_probs=get_emp_tip_probs(rf_states, tip_states, NaN, NaN)

#test=(hcat(tree.tlab,forbidden_aff, rf_states[tip_states], [rf_states[tip_probs[i].==1.0] for i in eachindex(tip_probs)]))


#rf_states[tip_states[30]]
#test=(hcat(tree.tlab, rf_states[tip_states], [rf_states[tip_probs[i].==1.0] for i in eachindex(tip_probs)]))
#
#show(stdout, "text/plain", hcat(eachindex(test), test[sortperm(test[:, 1]), :]))
#show(stdout, "text/plain", hcat(eachindex(test[:,1]), test))


#show(IOContext(io, :limit=>false), "text/plain", test)

#sum(sum.(tip_probs))/sum(sum.(full_uncertain_tip_probs))


#io = open((dir_name*"/"*"sim_info"*".txt"), "a") 
#writedlm(io,([file_name*"\t"*join(string.(rate_pars_sim).*"\t")]))                                                                                                                                                                                                                                                                                                                      
#close(io)        
#
#io = open((dir_name*"/"*"anc_states"*".txt"), "a") 
#writedlm(io,([file_name*"\t"*join(string.(anc_states).*"\t")]))                                                                                                                                                                                                                                                                                                                      
#close(io)        


if !isdir(dir_name*"/tip_states_sims")
    mkdir(dir_name*"/tip_states_sims")
    
end

##io = open((dir_name*"/tip_states/"*file_name*"tip_states"*".txt"), "a") #

#    writedlm(io,([join(string.(tip_states).*"\t")]))                                                                                                                                                                                                                                                                                                                      
#close(io)        



start_rates =[rand(rate_prior_dists[i]) for i in eachindex(move_types)]
start_clado_split_prob=[rand(clado_prior_dists[i]) for i in eachindex(clado_types)]

#join(round.(start_rates, sigdigits=3),"_" )

log_filename=(dir_name*"/"*file_name*"_log")




    realfun_mcmc(start_rates, clado_probs_start,
                 iters, iter_trims, write_interval,
                 anc_state_sampling, 
                 par_matrix, index_vec, move_types,rate_prior_dists,
                cladoPmat_unpar, clado_prior_dists, 
                 rf_states, 
                 tip_probs, 
                 tree, 
                 log_filename,
                 rate_tuning_par, clado_tuning_par,
                 rate_par_proposal_prob, 
                 prior_only,
                 false )




                 








##


burnin=0.5

#
#if !isdir(dir_name*"/post_pred")
#    mkdir(dir_name*"/post_pred")
#end
#
##RFBS_post_pred_log_filename=(dir_name*"/post_pred/"*par_chain_files[1]*"post_pred")
#
#
#RFBS_post_pred_log_filename=(dir_name*"/post_pred/"*"RFBS_"*file_name*"_post_pred")
#
##DEC posterior prediction            
#RFBS_log_mat = readdlm(log_filename)
#RFBS_par_chain=RFBS_log_mat[Int(round(size(RFBS_log_mat )[1]*burnin)):end,6:(6+length(move_types))]
#RFBS_model_par_names=RFBS_log_mat[1,6:(6+length(move_types))]
#
#
#RFBS_posterior_vec=[Float64.(RFBS_par_chain[i,:]) for i in 1:size(RFBS_par_chain)[1]]
#  
#
#RFBS_job_Bool_vec=[
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
#posterior_tip_pred_wrapper(post_pred_iters, post_pred_write_interval, 
#                            RFBS_posterior_vec, 
#                            tip_states, 
#                            tree_file,
#                           dir_name,
#                           RFBS_model_par_names,
#                           nbiomes,
#                           max_range,
#                           RFBS_post_pred_log_filename, false,
#                           RFBS_job_Bool_vec
#                           )      
#
#
#
#
#Q_zeros =zeros(length(rf_states), length(rf_states))
##file_name=string(rate_pars_sim)
#
#     
#
#prunelL_c, edge_probs_c, post_clado_probs_c = rfprune(rate_pars_sim     ,
#                     par_matrix       ,
#                     Q_zeros          ,
#                     index_vec        ,
#                     tip_probs        ,
#                     cladoPmat       ,
#                     brs              ,
#                     node_path        
#                              )
#priorlL_c = sum([log(pdf(prior_dists[i], rate_pars_sim[i])) for i in eachindex(rate_pars_sim)]) 
#lL_c= priorlL_c + prunelL_c  
##write files
#
#
#io = open((dir_name*"/"*"sim_info"*".txt"), "a") 
#writedlm(io,([file_name*"\t"*string(lL_c)*"\t"* string(priorlL_c) *"\t"* string(prunelL_c)*"\t"*join(string.(rate_pars_sim).*"\t")]))                                                                                                                                                                                                                                                                                                                      
#close(io)        
#
#io = open((dir_name*"/"*"rf_states"*".txt"), "a") 
#
#[writedlm(io,([i]))     for i in rf_states]                                                                                                                                                                                                                                                                                                                 
#close(io)        


#end 

#cat([estim_pars[i].-sim_pars[i] for i in eachindex(estim_pars)],dims=2)

#fund_adj_matrix
#forbidden_aff
#tip_probs_for=get_emp_tip_probs(rf_states, tip_states, NaN, forbidden_aff)
#
#tip_probs_for
#
#probs_comp=[sum(tip_probs_for[i].==tip_probs[i]) for i in eachindex(tip_probs)]
#probs_comp[probs_comp.<19]

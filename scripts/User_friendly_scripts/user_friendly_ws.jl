
cd(@__DIR__)

cd("../..")
pwd()

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


#include files 
include(joinpath(pwd(), "src", "CladoProbMatrix_fns.jl"))
include(joinpath(pwd(), "src", "make_rf_par_matrix_fns_switch.jl"))
include(joinpath(pwd(), "src", "Tree_Util_fns.jl"))
include(joinpath(pwd(), "src", "rf_sim_fns.jl"))
include(joinpath(pwd(), "src", "rf_rate_clado_mcmc_fns.jl"))
include(joinpath(pwd(), "src", "post_pred_fns.jl"))


biome_exc_treatment_ind = 2
biome_inc_treatment_ind = 2
adj_matrix             = true




data_dir        = joinpath("data", "emp", "viburnum_for_users")

reps=1
#start_truesim=true

biome_inc_treatment = ["cons" "bold" "none"][biome_inc_treatment_ind]
biome_exc_treatment = ["cons" "bold" "none"][biome_exc_treatment_ind]




# viburnum_dat=readdlm("viburnum_data_files/incf_viburnum_sorted_rf_states_3b.txt", Int64)
viburnum_dat = readdlm(joinpath("data", "emp", "viburnum", "included_affinity_treatments", biome_inc_treatment * "_incf_viburnum_sorted_rf_states_3b.txt"), Int64)
tree_file = joinpath(data_dir, "out.1.t163.f5.mask_fossil_states.mcc.tre") #pollen fossil taxa
 
tip_rf_ranges=[viburnum_dat[i,:] for i in 1:length(viburnum_dat[:,1])]

fn = joinpath(data_dir, "affs", biome_exc_treatment * "_excluded_" * biome_inc_treatment * "_included_combined_affinities.csv")
#readdlm("data/emp/viburnum_for_vignette/" * biome_exc_treatment * "_excluded_" * biome_inc_treatment     * "_included" * "_combined_affinities.csv", ",",String, header=true)


fund_adj_matrix=[0 1 0;  #biome 1 adjacencies
                 1 0 1 ;  #biome 2 adjacencies
                 0 1 0 ]  #biome3  adjacencies

tree, bts = read_tree(tree_file, order  = "cladewise",  branching_times = true, nexus_file_type=true) 


#RFBS arguements

data_fn               = fn
tree                  = tree
max_range             = 3
adj_matrix            = fund_adj_matrix
iters                 =1000000
iter_trims        =iters/10000
anc_state_sampling=iters/10000
write_interval           =100
post_pred_iters         =10000
post_pred_write_interval =1000
proposal_probs           = (rate_move = 10 , clado_move = 5 , rate_zero_switch = 2.0, clado_zero_switch = 0.0)
rate_tuning_par          =1.5
clado_tuning_par         =0.4

exponential_prior_rate   =0.5

allow_double_gains      = true
allow_double_losses     = false
allow_single_gains      = true
allow_single_losses     = true
allow_real_switches     = true
allow_realfun_switches  = false
ecological               =true
allopatric               =true


by_biome                = false
by_rf                   = true
by_gainloss             = true
by_doublesingle         = true




iters                 =1000
RFBS(
     data_fn               = fn
    ,tree                  = tree
    ,max_range             = 3
    ,adj_matrix            = fund_adj_matrix
    ,iters                 =1000
    ,iter_trims        =iters/100
    ,anc_state_sampling=iters/100
    ,write_interval           =100
    ,proposal_probs           = (rate_move = 10 , clado_move = 5 , rate_zero_switch = 2.0, clado_zero_switch = 0.0)
    ,rate_tuning_par          =1.5
    ,clado_tuning_par         =0.4
    ,exponential_prior_rate   =0.5
    ,allow_double_gains      = true
    ,allow_double_losses     = false
    ,allow_single_gains      = true
    ,allow_single_losses     = true
    ,allow_real_switches     = true
    ,allow_realfun_switches  = false
    ,ecological               =true
    ,allopatric               =true
    ,by_biome                = false
    ,by_rf                   = true
    ,by_gainloss             = true
    ,by_doublesingle         = true
    ,output_dirname           = joinpath(pwd(), "RFBS_user_example_output")
)
    


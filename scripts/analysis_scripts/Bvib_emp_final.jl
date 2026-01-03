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
include(pwd() *"/src/CladoProbMatrix_fns.jl")
include(pwd() *"/src/make_rf_par_matrix_fns_switch.jl")
include(pwd() *"/src/Tree_Util_fns.jl")
include(pwd() *"/src/rf_sim_fns.jl")
include(pwd() *"/src/rf_rate_clado_mcmc_fns.jl")
include(pwd() *"/src/post_pred_fns.jl")
include(pwd() *"/scripts/scratch_ws/RJ_update_mcmc_ws.jl")




##which are non zero Q matrix elements
allow_double_gains      = parse(Bool, ARGS[1])
allow_double_losses     = parse(Bool, ARGS[2])
allow_single_gains      = parse(Bool, ARGS[3])
allow_single_losses     = parse(Bool, ARGS[4])
by_biome                = parse(Bool, ARGS[5])
by_rf                   = parse(Bool, ARGS[6])
by_gainloss             = parse(Bool, ARGS[7])
by_doublesingle         = parse(Bool, ARGS[8])
DEC                     = parse(Bool, ARGS[9])
include_fund            = parse(Bool, ARGS[10])
use_forbidden_fund_vec  = parse(Bool, ARGS[11])
use_adj_matrix          = parse(Bool, ARGS[12])
RJ_ana                  = parse(Bool, ARGS[13]) #if true use RJ update for clado probs, if false use MCMC update
RJ_clado                = parse(Bool, ARGS[14]) #if true use RJ update for clado probs, if false use MCMC update

biome_exc_treatment_ind = parse(Int, ARGS[15])
biome_inc_treatment_ind = parse(Int, ARGS[16])

run_num                 = parse(Int, ARGS[17])



allow_real_switches    = true
allow_realfun_switches = false


data_dir        = "data/emp/viburnum/"
exclude_aff_dir = data_dir*"excluded_aff/"
include_aff_dir = data_dir*"included_aff/"

reps=1
#start_truesim=true
Fossil=true

biome_inc_treatment = ["cons" "bold"][biome_inc_treatment_ind]

if include_fund


   # viburnum_dat=readdlm("viburnum_data_files/incf_viburnum_sorted_rf_states_3b.txt", Int64)
   viburnum_dat=readdlm(include_aff_dir  * biome_inc_treatment * "_incf_viburnum_sorted_rf_states_3b.txt", Int64)

    if Fossil
        tree_file = data_dir * "out.1.t163.f5.mask_fossil_states.mcc.tre" #pollen fossil taxa
    else

        tree_file = data_dir * "viburnum_sorted.tre" #extant only
    end
    tip_rf_ranges=[viburnum_dat[i,:] for i in 1:length(viburnum_dat[:,1])]


else

    viburnum_dat = readdlm( data_dir * "viburnum_sorted_rf_states_3b.txt", Int64)
    viburnum_dat=readdlm(include_aff_dir  * biome_inc_treatment * "_incf_viburnum_sorted_rf_states_3b.txt", Int64)
    viburnum_dat[(viburnum_dat.==1)] .= 0

    if Fossil
        tree_file = data_dir*"out.1.t163.f5.mask_fossil_states.mcc.tre" #pollen fossil taxa
    else

        tree_file = data_dir*"viburnum_sorted.tre" #extant only
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
    
         forbidden_aff = Int.(readdlm(exclude_aff_dir* biome_exc_treatment*".csv", ',', header=true)[1])
    
 
    end

else
    forbidden_aff=NaN
end


if RJ_ana
    rate_zero_switch=2.0
else
    rate_zero_switch=0.0
end
if RJ_clado
    clado_zero_switch=2.0
else
    clado_zero_switch=0.0
end
iters                 =1000000


iter_trims        =iters/10000
anc_state_sampling=iters/10000
write_interval           =100
post_pred_iters         =10000
post_pred_write_interval =1000
proposal_probs           = (rate_move = 10 , clado_move = 5 , rate_zero_switch = rate_zero_switch, clado_zero_switch = clado_zero_switch)
rate_tuning_par          =1.5
clado_tuning_par         =0.4
ecological               =true
allopatric               =true

if DEC

    ecological=true
    allopatric=false

end

Prior=0.5
rate_prior_vec  = [Exponential(Prior)]

clado_types=["sub_split"]


prior_label="Exp" * replace(replace(string(Prior),"."=> "p"), "-" => "neg")

prior_only=true

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



dir_name = "outfiles/emp/viburnum/resub/" * "Bvib_"* string(nbiomes_g[1])*"nB_"*prior_label*"_"* string(iters)


if prior_only==true
    dir_name=dir_name*"_PO"
end

if Fossil==true
    dir_name=dir_name*"_Foss"
end


if use_adj_matrix==true
    dir_name=dir_name*"_admat"
end

if include_fund==true
    dir_name=dir_name*"_incf"
    dir_name=dir_name * "_" * biome_inc_treatment
    
else
    dir_name=dir_name*"_noincf"
    dir_name=dir_name * "_" * biome_inc_treatment
end


if use_forbidden_fund_vec==true
    dir_name=dir_name*"_excf"
    dir_name=dir_name * "_" * biome_exc_treatment

else
    dir_name=dir_name*"_noexcf"
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

if allow_real_switches==true
    dir_name=dir_name*"_2sw"
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


if Fossil==true
    dir_name=dir_name*"_Foss"
end

if RJ_ana==true
    dir_name=dir_name*"_RJa"
end
if RJ_clado==true
    dir_name=dir_name*"_RJc"
end

dir_name = replace(dir_name,r" "  => s"_")


if !isdir(dir_name)
    mkdir(dir_name)
    
end


estim_pars=[]
Q_set=[]
start_sets=[]
out_loglik_set=[]

NaN_run=[]



#
print("rep")
print(" ")
#print(i)
print(" ")

nbiomes=nbiomes_g[1]

max_range=max_range_g[1]

allow_realfun_switches = false

# take details to generater Q vectors to map rates to right Q matrix, the Q matrix generated by this function is ignored (just using the simulation function, messy)
par_matrix, rf_states, move_types, index_vec, move_matrix = make_rf_par_matrix(nbiomes,max_range,
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


rf2stDict = Dict(copy(rf_states).=>1:length(rf_states))
st2rfDict = Dict( (1:length(rf_states)).=>copy(rf_states))

tip_states=[rf2stDict[tip_rf_ranges[i]] for i in eachindex(tip_rf_ranges)]

cladoPmat_unpar= makeclado_Pmat_equal_prob(rf_states,
ecological,
allopatric)


file_name="Bvib__"*string(run_num)*"_"


cladoPmat_sim, clado_probs_start, split_index_vec, sub_index_vec, String_Clado_Mats, clado_prior_dists =make_clado_Pmat(rf_states, ecological, allopatric )


rate_prior_dists=make_Prior(rate_prior_vec, move_types)


tree, bts = read_tree(tree_file, order  = "cladewise",  branching_times = true, nexus_file_type=true) 
ntips=length(tree.tlab)
ed=tree.ed
edges = cat(ed, [2*ntips ntips + 1], dims = 1)
triads=maketriads(edges)
node_path  = get_trav_path_4prun(triads )
br = branching_times(tree)
brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
 

tip_probs = get_emp_tip_probs(rf_states, tip_states, fund_adj_matrix, forbidden_aff)



tree.tlab
rf_states
full_uncertain_tip_probs=get_emp_tip_probs(rf_states, tip_states, NaN, NaN)


if !isdir(dir_name*"/tip_states_sims")
    mkdir(dir_name*"/tip_states_sims")
    
end



start_rates =[rand(rate_prior_dists[i]) for i in eachindex(move_types)]
start_clado_split_prob=[rand(clado_prior_dists[i]) for i in eachindex(clado_types)]


log_filename=(dir_name*"/"*file_name*"_log")




    realfun_mcmc(start_rates,                                              
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
                 tip_probs, 
                 tree, 
                 log_filename,
                 rate_tuning_par, 
                 clado_tuning_par,
                 proposal_probs, 
                 prior_only,
                 false )




##functions for specifying and generating a realfun Q matrix


function convert_RF_2_DEC_states(rf_state_vec, rf_states, DEC_states, fun2non=true)

    DEC_state_vec=copy(rf_state_vec)

    if fun2non
        RF_real_aff=[rf_states[i].==2 for i in eachindex(rf_states)]
    else
        RF_real_aff=[rf_states[i].>0 for i in eachindex(rf_states)]
    end

    DEC_real_aff=[DEC_states[i].==2 for i in eachindex(DEC_states)]

    #create a n=length of RFBS state space mapping RFBS states to new DEC states
    rf2dec =   [findall(x->x==1,[RF_real_aff[i]].==DEC_real_aff)[1] for i in eachindex(RF_real_aff)]


    #Idfentify the one zero in the state vec, this is the separator between the node and cladogenetic ancestral states
     empty_state_ind=findall(x->x==0,rf_state_vec)

    for i in eachindex(rf_state_vec)

        #print(i)
        if !in(i,  empty_state_ind)

            DEC_state_vec[i]=rf2dec[rf_state_vec[i]]

        end
    end


    return DEC_state_vec

end


function make_DEC_states(nbiomes::Int64,max_real_biome::Int64)

    permutation=collect(combinations(reduce(vcat,[[0,2] for i in 1:nbiomes])))
  
    nbiome_states=unique(permutation[(length.(permutation).==nbiomes)] )
  
    all_realfun_biome_states=nbiome_states[[2 in i for i in nbiome_states]]
  
    #realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==1) <= 0) for i in 1:length(all_realfun_biome_states) ]]
  
    realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==2) <= max_real_biome) for i in 1:length(all_realfun_biome_states) ]]


    realfun_biome_states =sort(realfun_biome_states)
  
    return(realfun_biome_states)
end


#generate vector of possible biome ranges a lineage can occupy (ex, for two biomes with one fundamental and realized affinity you have [1,2] )
function make_rf_states(nbiomes::Int64,max_real_biome::Int64)

    permutation=collect(combinations(reduce(vcat,[[0,1,2] for i in 1:nbiomes])))
  
    nbiome_states=unique(permutation[(length.(permutation).==nbiomes)] )
  
    all_realfun_biome_states=nbiome_states[[2 in i for i in nbiome_states]]
  
    #realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==1) <= 0) for i in 1:length(all_realfun_biome_states) ]]
  
    realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==2) <= max_real_biome) for i in 1:length(all_realfun_biome_states) ]]


    realfun_biome_states =sort(realfun_biome_states)
  
    return(realfun_biome_states)
end


#generate vector of possible biome ranges a lineage can occupy (ex, for two biomes with one fundamental and realized affinity you have [1,2] )
function make_rf_states(nbiomes::Int64,max_real_biome::Int64, no_fun::Bool=false)

    permutation=collect(combinations(reduce(vcat,[[0,1,2] for i in 1:nbiomes])))
  
    nbiome_states=unique(permutation[(length.(permutation).==nbiomes)] )
  
    all_realfun_biome_states=nbiome_states[[2 in i for i in nbiome_states]]
  
    if no_fun==true
        all_realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==1) <= 0) for i in 1:length(all_realfun_biome_states) ]]
    end

    realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==2) <= max_real_biome) for i in 1:length(all_realfun_biome_states) ]]


    realfun_biome_states =sort(realfun_biome_states)
  
    return(realfun_biome_states)
end


#generate tip probs based on tip states etiehr with full certainty or a sampled percentage of uncertainty for each uncertain state, equal tip probability is given to, the realized affinity with no fundamental, the actual simulated tip state and a set percentage of "droppable" states that are neither the simulated state nor the realized affinity without fundamental 
function  get_tip_probs(rf_states, tip_states ,uncertain_tips, uncertain_percent, down_sample_state_group_percent, state_group_percent,  drop_real)
              
    #find index for real biome affinities in each range
    
    if uncertain_tips==true

        real_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
        real_groups=unique(real_ind)
        state_groups=[findall(item -> item == i, real_ind) for i in real_groups]
    
   
        tip_probs=Vector{Vector{Float64}}(undef, length(tip_states))
        uncertain_tips_vec=sort(sample(eachindex(tip_states),Integer(round(length(tip_states)*uncertain_percent)), replace=false))
        state_group_tips_vec=sort(sample(uncertain_tips_vec,Integer(round(length(tip_states)*down_sample_state_group_percent)), replace=false))


        for i in eachindex(tip_states)
            #push!(tip_probs, fill(0.0,length(rf_states)) )
            tip_probs[i]=fill(0.0,length(rf_states))
            #initializes vectors of probabilities and fills o0nes in for each state in the tips state group, equal probs for all putative states, how to paramaterize this?
       
            if in(i, uncertain_tips_vec) 


                state_group= state_groups[in.(tip_states[i],state_groups)][1]


                #essential_states = append!(state_group[1 .∉ rf_states[state_group]],tip_states[i])
                
                if drop_real

                    essential_states = [tip_states[i]]
                    else

                    essential_states = append!(state_group[1 .∉ rf_states[state_group]],tip_states[i])
                 end

                droppable_states=setdiff(state_group, essential_states)

                Nstates_kept=Integer(ceil(length(droppable_states)*state_group_percent))

                droppable_subsample=sample(droppable_states, Nstates_kept, replace=false)


                state_group_samp=unique(append!(essential_states, droppable_subsample))


               if  in( i,  state_group_tips_vec) 

                tip_probs[i]= Float64.(1.0*(in.((eachindex(rf_states)),[state_group_samp]))) 

               else

                tip_probs[i]= Float64.(1.0*(in.((eachindex(rf_states)),[state_group]))) 



               end

                #sample from the state group (while keeping the realized only state)
                #state_group_samp=unique(append!(state_group[sort(sample(eachindex(state_group),Integer(ceil(length(state_group)*state_group_percent)), replace=false))], state_group[1 .∉ rf_states[state_group]]))
                
            else
                tip_probs[i]=  Float64.(1.0*(in.((eachindex(rf_states)), tip_states[i])) )
            end
        end

    else
    
        tip_probs=[  Float64.(1.0*(in.((eachindex(rf_states)), tip_states[i]))) for i in eachindex(tip_states)]
    end
    
    
    #if uncertain_tips==true
    ##initializes vectors of probabilities and fills o0nes in for each state in the tips state group, equal probs for all putative states, how to paramaterize this?
    #    tip_probs=[ 1.0*(in.((eachindex(rf_states)),state_groups[in.(tip_states[i],state_groups)])) for i in 1:length(tip_states)]
    #else
    #    tip_probs=[ 1.0*(in.((eachindex(rf_states)), tip_states[i])) for i in 1:length(tip_states)]
    #end
    
    #return(tip_probs::Vector{Vector{Float64}},state_groups,  state_group_tips_vec)

    return(tip_probs::Vector{Vector{Float64}})

end


function drop_fund_adj_states(obs_state, #observed state
                             state_group, #group of possible states given observed state
                             fund_adj_matrix, #adjacency matrix, for each state row, shows which state columns are allowed or not
                             rf_states #state space
)


    fund_allowed_by_biome=fund_adj_matrix[rf_states[obs_state].>1,:]
    fund_allowed=sum(fund_allowed_by_biome,dims=1).>0

    new_state_groups=[]

    for group_state in state_group[obs_state.∉ state_group]

        #group_state=state_group[obs_state.∉ state_group][1]
        #group_state=14

        #save only indices that are not realized but allowed as fund, we want to use this to check which state have the fundamental affinity from the state group
        bad_fund=findall(item -> fund_allowed[item]==0 ,1:length(rf_states[1]))

        #bridge fund are the fundamental affinities needed for a fundamental affinity in the bad fund to be possible (ex for tropical to have a fundamental affinity for cold it must also have an affinity for warm temperate)
       

        #rf_states[group_state][bad_fund]==fund_allowed[bad_fund]
    
        #get all possible fun affinities for this state
        #all_fund=setdiff(1:nbiomes, real_ind[group_state])
       # all_fund=rf_states[real_ind[state]]==
    

       #include state if none of the bad funds are present (they are all zeros) or if they are ALL present 
        if !in(1, rf_states[group_state][bad_fund]) || (all(rf_states[group_state][bad_fund].>0)&& all(rf_states[group_state][vec(fund_allowed)].>0)) 
            #tip_prob_mat[state,group_state]=1
            push!(new_state_groups,group_state)
        
        end
    
       # tip_prob_mat[state,state]=1
    
    end 

    push!(new_state_groups,obs_state)

    return(new_state_groups)

end


function drop_forbid_fund_states(obs_state, #observed state
                               state_group, #group of possible states given observed state
                               sp_forbid_aff, #adjacency matrix, for each state row, shows which state columns are allowed or not
                               rf_states #state space
)


    new_state_groups=[]
    for group_state in state_group[obs_state.∉ state_group]
        #save only indices that are not realized but allowed as fund, we want to use this to check which state have the fundamental affinity from the state group
        bad_fund=(1:length(rf_states[1])).==sp_forbid_aff
    
        #rf_states[group_state][bad_fund]==fund_allowed[bad_fund]
    
        #get all possible fun affinities for this state
        #all_fund=setdiff(1:nbiomes, real_ind[group_state])
       # all_fund=rf_states[real_ind[state]]==
    
        if !in(1, rf_states[group_state][bad_fund])
            #tip_prob_mat[state,group_state]=1
            push!(new_state_groups,group_state)
        
        end
    
       # tip_prob_mat[state,state]=1
    
    end 

    push!(new_state_groups,obs_state)


    return(new_state_groups)

end


#generate tip probs based on tip states and possible unobserved fund affinities, can include adjacency mappings to see which fund aff are prohibited with certain realized aff


function   get_emp_tip_probs(rf_states, tip_states, fund_adj_matrix)

    real_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
    real_groups=unique(real_ind)
    state_groups=[findall(item -> item == i, real_ind) for i in real_groups]

    tip_prob_mat=fill(0,length(rf_states), length(rf_states))


    tip_probs=Vector{Vector{Float64}}(undef, length(tip_states))


    for i in eachindex(tip_states)
        #push!(tip_probs, fill(0.0,length(rf_states)) )
        tip_probs[i]=fill(0.0,length(rf_states))

        fund_allowed_mat=fund_adj_matrix[real_ind[tip_states[i]],:]

        fund_allowed=[in(1, fund_allowed_mat[:,row]) for row in 1:length(fund_allowed_mat[1,:])]

        #initializes vectors of probabilities and fills o0nes in for each state in the tips state group, equal probs for all putative states, how to paramaterize this?



        state_group=state_groups[in.(tip_states[i],state_groups)][1]



        if(1 in rf_states[tip_states[i]])

            state_group=state_group[[rf_states[tip_states[i]].==1].==[rf_states[st].==1 for st in state_group]]


        end
        tip_probs[i]= Float64.(1.0*(in.((eachindex(rf_states)),[state_group]))) 

        
      

        #essential_states = append!(state_group[1 .∉ rf_states[state_group]],tip_states[i])


        for group_state in state_group[tip_states[i].∉ state_group]

            #save only indices that are not realized but allowed as fund, we want to use this to check which state have the fundamental affinity from the state group
            bad_fund=findall(item -> fund_allowed[item]==0 ,1:length(rf_states[1]))
        

            rf_states[group_state][bad_fund]==fund_allowed[bad_fund]
        
            #get all possible fun affinities for this state
            #all_fund=setdiff(1:nbiomes, real_ind[group_state])
           # all_fund=rf_states[real_ind[state]]==
        

            if in(1, rf_states[group_state][bad_fund])
                #tip_prob_mat[state,group_state]=1
                tip_probs[i][group_state]=0.0
            
            end
        
           # tip_prob_mat[state,state]=1
        
        end 


        #sample from the state group (while keeping the realized only state)
        #state_group_samp=unique(append!(state_group[sort(sample(eachindex(state_group),Integer(ceil(length(state_group)*state_group_percent)), replace=false))], state_group[1 .∉ rf_states[state_group]]))

    end

    return(tip_probs)

end

function   get_emp_tip_probs(rf_states, tip_states, fund_adj_matrix=NaN, forbidden_fund=NaN)


    real_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
    real_groups=unique(real_ind)
    state_groups=[findall(item -> item == i, real_ind) for i in real_groups]

    #tip_prob_mat=fill(0,length(rf_states), length(rf_states))


    tip_probs=Vector{Vector{Float64}}(undef, length(tip_states))


    for i in eachindex(tip_states)
        #push!(tip_probs, fill(0.0,length(rf_states)) )
        print(i)
        tip_probs[i]=fill(0.0,length(rf_states))

       # full_tip_state_group=state_groups[in.(tip_states[i],state_groups)][1]

        tip_state_group=state_groups[in.(tip_states[i],state_groups)][1]

        #tip_state_group= state_groups[tip_states[i]]
        if(1 in rf_states[tip_states[i]])

            tip_state_group=tip_state_group[[rf_states[tip_states[i]].==1].==[rf_states[st].==1 for st in tip_state_group]]


        end



        if fund_adj_matrix!==NaN
            tip_state_group=drop_fund_adj_states(tip_states[i], #observed state
                                                tip_state_group, #group of possible states given observed state
                                                fund_adj_matrix, #adjacency matrix, for each state row, shows which state columns are allowed or not
                                                rf_states #state space
            )

        end


        if forbidden_fund!==NaN


            if(size(forbidden_aff)[2]>1)

                for col in 1:length(size(forbidden_aff))

                    tip_state_group=drop_forbid_fund_states(tip_states[i], #observed state
                    tip_state_group, #group of possible states given observed state
                    forbidden_fund[i, col], #adjacency matrix, for each state row, shows which state columns are allowed or not
                    rf_states #state space
                    )

                end

            else 
                tip_state_group=drop_forbid_fund_states(tip_states[i], #observed state
                                                    tip_state_group, #group of possible states given observed state
                                                    forbidden_fund[i], #adjacency matrix, for each state row, shows which state columns are allowed or not
                                                    rf_states #state space
                                                    )

            end
        end
        
        
       #fund_aff_bit=rf_states[tip_states[i]].==1
        
       rf_states[tip_state_group]

        #state_group_samp=unique(append!(state_group[sort(sample(eachindex(state_group),Integer(ceil(length(state_group)*state_group_percent)), replace=false))], state_group[1 .∉ rf_states[state_group]]))
        tip_probs[i][tip_state_group].=1.0
        tip_probs[i][tip_states[i]]=1.0



    end

    return(tip_probs)

end

#generate uncertain tip probs based on vector of tips to drop certainty on
function  make_uncertain_tips(rf_states,tip_states, uncertain_tips_vec)

    real_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
    real_groups=unique(real_ind)
    state_groups=[findall(item -> item == i, real_ind) for i in real_groups]


    for i in eachindex(tip_states)
        #push!(tip_probs, fill(0.0,length(rf_states)) )
        tip_probs[i]=fill(0.0,length(rf_states))
        #initializes vectors of probabilities and fills o0nes in for each state in the tips state group, equal probs for all putative states, how to paramaterize this?

        if in(i, uncertain_tips_vec) 
            tip_probs[i]= Float64.(1.0*(in.((eachindex(rf_states)),state_groups[in.(tip_states[i],state_groups)]))) 
        else
            tip_probs[i]=  Float64.(1.0*(in.((eachindex(rf_states)), tip_states[i])) )
        end
    end

    return(tip_probs)

end

#generate uncertain tip probs based on vector of tips to drop certainty on


#sort(sample(eachindex(tip_states),Integer(round(length(tip_states)*uncertain_percent)), replace=false))

# extracts  rates from matrix of parameter strings to match up integer move types with an empty integer rate matrix (should result in a corHMM/phytools matrix )
function map_intpars2emptymat(par_matrix::Matrix{String},
                              par_matrix_numeric:: Matrix{Int64})

    moves_types=unique(par_matrix)[2:end]




    par_matrix_numeric=fill(0,length(par_matrix[1,:]), length(par_matrix[:,1]))

    for i in eachindex(moves_types)

        par_matrix_numeric[par_matrix.===moves_types[i]] .= i

    end

    return(par_matrix_numeric)

end


#generate an integer Q matrix where free pars are scored as integers, allow types of moves and specify types of free parameters 
function make_rf_par_matrix(nbiomes::Int64,
                            max_real_biome::Int64,
                            allow_double_gains::Bool  =true,
                            allow_double_losses::Bool =false,
                            allow_single_gains::Bool  =true,
                            allow_single_losses::Bool =true,
                            by_biome           ::Bool =true,
                            by_rf              ::Bool =true,
                            by_gainloss        ::Bool =true,
                            by_doublesingle    ::Bool =true,
                            DEC                ::Bool=false
                            )


    if(DEC==true)
        #if DEC model then no fundamental affinities are modeled
        no_fun=true

        biome_range_states=make_rf_states(nbiomes, max_real_biome, no_fun)
        #if DEC all double moves must be allowed
        allow_double_gains=true
        allow_double_losses=true

        #free pars for double vs single moves and real vs fundamental moves arent allowed 
        by_rf =false
        by_doublesingle=false
    else    
        biome_range_states=make_rf_states(nbiomes, max_real_biome)
    end
    #biome_range_states=make_rf_states(3, 2)

    #final matrix that will have zeros at non rates and 1 + for real rates
    par_matrix_int_empty=fill(0,length(biome_range_states), length(biome_range_states))

    par_vec=[]

    #create array of strings to put rates in (easier than trying to do it directly to numbers )
    par_matrix_string= fill("-", length(biome_range_states), length(biome_range_states))
    
    #iterate throug each element in Q matrix if and based on moves allowed 
    #different changing invalid moves from "-" to a "MV" with more details 
    #added if more than one free par is specified

    for from in eachindex(biome_range_states) for to in eachindex(biome_range_states) 

        biome_diff = biome_range_states[to]-biome_range_states[from]
  
        #find shifts that equal have only single gains or loses
        #single_shifts=(biome_diff .== 1 .||  biome_diff .== -1)
        single_gains=(biome_diff .== 1 )
        single_losses=(biome_diff .== -1 )
  
        #check if it is a double gain 
        double_gains=(biome_diff .== 2 )
        double_losses=(biome_diff .== -2 )
  
  
        no_shifts=(biome_diff .== 0) 
  
  
        #check to make sure there is only a shift in a single biome, those are the only possible valid moves
        if (sum(no_shifts) == (nbiomes-1)) 
  
          
          #check if a single rf gain and then fill other details
          if (sum(single_gains) == (1)) && allow_single_gains==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[to][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
              
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_g" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_s"
              end 
  
              push!(par_vec, [from, to])
  
          end
     
          #check if a single rf biome loss and then fill other details
          if (sum(single_losses) == (1)) && allow_single_losses==true
    
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[from][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
              
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_l" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_s"
              end 
    
              push!(par_vec, [from, to])
  
          end
  
          #check if a double rf biome gains and then fill other details
          if (sum(double_gains) == (1)) && allow_double_gains==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[to][biome_diff.!=0][1])
              end
          
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
  
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_g" 
              end
          
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_d"
              end 
  
              push!(par_vec, [from, to])
  
          end
  
          #check if a double rf biome loss and then fill other details
          if (sum(double_losses) == (1)) && allow_double_losses==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[from][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
  
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_l" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_d"
              end 
  
              push!(par_vec, [from, to])
  
          end
  
        end 
  
    end end 
  
    par_matrix_int = map_intpars2emptymat(par_matrix_string, par_matrix_int_empty)

    moves_types=unique(par_matrix_string)[2:end]

    return(par_matrix_int,  biome_range_states, moves_types, par_vec, par_matrix_string)

end


#generate an integer Q matrix where free pars are scored as integers, allow types of moves and specify types of free parameters 
function make_rf_par_matrix(nbiomes::Int64,
                            max_real_biome::Int64,
                            allow_double_gains::Bool  =true,
                            allow_double_losses::Bool =false,
                            allow_single_gains::Bool  =true,
                            allow_single_losses::Bool =true,
                            by_biome           ::Bool =true,
                            by_rf              ::Bool =true,
                            by_gainloss        ::Bool =true,
                            by_doublesingle    ::Bool =true
                            )

    biome_range_states=make_rf_states(nbiomes, max_real_biome)

    #biome_range_states=make_rf_states(3, 2)

    #final matrix that will have zeros at non rates and 1 + for real rates
    par_matrix_int_empty=fill(0,length(biome_range_states), length(biome_range_states))

    par_vec=[]

    #create array of strings to put rates in (easier than trying to do it directly to numbers )
    par_matrix_string= fill("-", length(biome_range_states), length(biome_range_states))
    
    #iterate throug each element in Q matrix if and based on moves allowed 
    #different changing invalid moves from "-" to a "MV" with more details 
    #added if more than one free par is specified

    for from in eachindex(biome_range_states) for to in eachindex(biome_range_states) 

        biome_diff = biome_range_states[to]-biome_range_states[from]
  
        #find shifts that equal have only single gains or loses
        #single_shifts=(biome_diff .== 1 .||  biome_diff .== -1)
        single_gains=(biome_diff .== 1 )
        single_losses=(biome_diff .== -1 )
  
        #check if it is a double gain 
        double_gains=(biome_diff .== 2 )
        double_losses=(biome_diff .== -2 )
  
  
        no_shifts=(biome_diff .== 0) 
  
  
        #check to make sure there is only a shift in a single biome, those are the only possible valid moves
        if (sum(no_shifts) == (nbiomes-1)) 
  
          
          #check if a single rf gain and then fill other details
          if (sum(single_gains) == (1)) && allow_single_gains==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[to][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
              
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_g" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_s"
              end 
  
              push!(par_vec, [from, to])
  
          end
     
          #check if a single rf biome loss and then fill other details
          if (sum(single_losses) == (1)) && allow_single_losses==true
    
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[from][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
              
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_l" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_s"
              end 
    
              push!(par_vec, [from, to])
  
          end
  
          #check if a double rf biome gains and then fill other details
          if (sum(double_gains) == (1)) && allow_double_gains==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[to][biome_diff.!=0][1])
              end
          
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
  
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_g" 
              end
          
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_d"
              end 
  
              push!(par_vec, [from, to])
  
          end
  
          #check if a double rf biome loss and then fill other details
          if (sum(double_losses) == (1)) && allow_double_losses==true
  
              par_matrix_string[from, to]=""
  
              if by_rf==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_rf" *  string(biome_range_states[from][biome_diff.!=0][1])
              end
  
              if by_biome==true
                  par_matrix_string[from, to]= par_matrix_string[from, to]* "_b"* string(eachindex(biome_diff)[biome_diff.!=0][1])  
              end
  
              if by_gainloss==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_l" 
              end
  
              if by_doublesingle==true
                  par_matrix_string[from, to]= par_matrix_string[from, to] * "_d"
              end 
  
              push!(par_vec, [from, to])
  
          end
  
        end 
  
    end end 
  
    par_matrix_int = map_intpars2emptymat(par_matrix_string, par_matrix_int_empty)

    moves_types=unique(par_matrix_string)[2:end]

    return(par_matrix_int,  biome_range_states, moves_types, par_vec, par_matrix_string)

end


function fill_emptyQmat(Q_mat_zeros  ::Matrix{Float64} ,
                        rate_pars     :: Vector{Float64},
                        par_vec       ::Vector{Any}     ,
                        par_matrix_int::Matrix{Int64}    
    )

    Q_mat_empty = deepcopy(Q_mat_zeros)         
    
    for i in par_vec
  
        Q_mat_empty[i[1],i[2]]=rate_pars[par_matrix_int[i[1],i[2]]]
  
    end


    for state in  1:length(Q_mat_empty[1,:])

        Q_mat_empty[state,state] = -sum(Q_mat_empty[state, 1:end ] )
    
    end
    
    Q_mat_full = Q_mat_empty

  return( Q_mat_full)
end











##generate vector of possible biome ranges a lineage can occupy (ex, for two biomes with one fundamental and realized affinity you have [1,2] )
#function make_rf_states(nbiomes::Int64,max_real_biome::Int64)
#
#    permutation=collect(combinations(reduce(vcat,[[0,1,2] for i in 1:nbiomes])))
#  
#    nbiome_states=unique(permutation[(length.(permutation).==nbiomes)] )
#  
#    all_realfun_biome_states=nbiome_states[[2 in i for i in nbiome_states]]
#  
#    realfun_biome_states=all_realfun_biome_states[[ (sum(all_realfun_biome_states[i].==2) <= max_real_biome) for i in 1:length(all_realfun_biome_states) ]]
#  
#    realfun_biome_states =sort(realfun_biome_states)
#
#    real_ind=[findall(item -> item == 2, realfun_biome_states[i]) for i in eachindex(realfun_biome_states)]
#    real_groups=unique(real_ind)
#    state_groups=[findall(item -> item == i, real_ind) for i in real_groups]
#    
#    tip_prob_mat=fill(0,length(realfun_biome_states), length(realfun_biome_states))
#    
#    filtered_realfun_biome_states=[]
#
#    for state in eachindex(realfun_biome_states)
#    
#    fund_allowed_mat = fund_adj_matrix[ real_ind[state] , : ]
#    
#       fund_allowed = [in( 1,fund_allowed_mat[ :,row ] ) for row in 1:length( fund_allowed_mat[ 1,:]  ) ]
#    
#        for group_state in state_groups[in.(state,state_groups)][1]
#    
#            #save only indices that are not realized and allowed as fund, we want to use this to check which state have the fundamental affinity from the state group
#            fund_ind=findall(item -> !in(item ,real_ind[state])&fund_allowed[item]==1, 1:length(realfun_biome_states[1]))
#    
#            realfun_biome_states[group_state][fund_ind]==fund_allowed[fund_ind]
#            
#            
#    
#            if realfun_biome_states[group_state][fund_ind]==fund_allowed[fund_ind]
#    
#        
#               push!(filtered_realfun_biome_states, realfun_biome_states[group_state])
#    
#            elseif sum(realfun_biome_states[group_state][fund_allowed])==0
#                tip_prob_mat[state,group_state]=1
#            end
#    
#            #tip_prob_mat[state,state]=1
#    
#
#
#        end 
#    end
#    
#    [sum(tip_prob_mat[:,col] for col in 1:length(rf_states))]
#
#    return(realfun_biome_states)
#end
#
##fund_cartind=[findall(item -> item == 1, fund_adj_matrix[real_ind[i],:])  for i in eachindex(real_ind)]
#
#
##fund_ind=[[fund_cartind[i][coord][2] for coord in eachindex(fund_cartind[i])] for i in eachindex(fund_cartind) ]
#
#
#
#
#
##fun_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
#
#
#real_ind=[findall(item -> item == 2, rf_states[i]) for i in eachindex(rf_states)]
#real_groups=unique(real_ind)
#state_groups=[findall(item -> item == i, real_ind) for i in real_groups]
#
#tip_prob_mat=fill(0,length(rf_states), length(rf_states))
#
#inclusive=true
#for state in eachindex(rf_states)
#
#    fund_allowed_mat=fund_adj_matrix[real_ind[state],:]
#
#   fund_allowed=[in(1, fund_allowed_mat[:,row]) for row in 1:length(fund_allowed_mat[1,:])]
#
#    for group_state in state_groups[in.(state,state_groups)][1]
#
#        #save only indices that are not realized but allowed as fund, we want to use this to check which state have the fundamental affinity from the state group
#        fund_ind=findall(item -> !in(item ,real_ind[state])&fund_allowed[item]==1, 1:length(rf_states[1]))
#
#        rf_states[group_state][fund_ind]==fund_allowed[fund_ind]
#        
#        rf_states[group_state][fund_ind]==fund_allowed[fund_ind]
#
#        all_fund=setdiff(1:nbiomes, real_ind[group_state])
#       # all_fund=rf_states[real_ind[state]]==
#
#        
#        if rf_states[group_state][all_fund]==fund_allowed[all_fund]
#
#    
#            tip_prob_mat[state,group_state]=1
#
#        elseif rf_states[group_state][fund_ind]==fund_allowed[fund_ind] && inclusive
#            tip_prob_mat[state,group_state]=1
#
#        elseif sum(rf_states[group_state][all_fund])==0
#            tip_prob_mat[state,group_state]=1
#        end
#
#       # tip_prob_mat[state,state]=1
#
#    end 
#end
#
#
#in(0, tip_prob_mat[:,col])
#
#filtered_rf_states=rf_states[[ sum(tip_prob_mat[:,col])>0 for col in 1:length(rf_states)]]
#
#
#fund_ind=findall(item -> in(item ,real_ind[state]), 1:3)
#
#rf_states[group_state][fund_ind]==fund_allowed[fund_ind]
#

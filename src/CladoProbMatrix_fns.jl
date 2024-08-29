
function get_other_child_range(p_range::Array{Int64, 1}, c1::Vector{Vector{Int64}}, retain_fun::Bool=true)

  #generates second child state for split ranges (where neither fully inherits the parent state)
  if retain_fun
      post_clado=1
  else
      post_clado=0
  end 

  c2 =deepcopy(c1)

  for i in eachindex(c1) for j in eachindex(c1[i])

    #print(i)
    #print(j)
  
      if c1[i][j]==2
      
        c2[i][j]=post_clado
      
      elseif  c1[i][j]==post_clado && p_range[j]!=post_clado
      
        c2[i][j]=2
      
      end 

     
  end  end

  return(c2)
end 



function get_all_clados(p_range::Array{Int64, 1}, realfun_states::Array{Array{Int64, 1}, 1})

  #generates second child state for split ranges (where neither fully inherits the parent state)


  #1[1] p_range
  #
  #ealfun_states[18][p_range.==0].== p_range[p_range.==0]
  #ealfun_states[18][p_range.>0].== p_range[p_range.>0]
  #1[1][p_range.>0].== p_range[p_range.>0]


  c1 = []
  c2 = []

  for i in eachindex(realfun_states) 

    for j in eachindex(realfun_states) 

      #print(i)
      #print(j)
 
      
     #c1_nonaff_equal = realfun_states[i][p_range.==0].== p_range[p_range.==0]
     #c2_nonaff_equal = realfun_states[j][p_range.==0].== p_range[p_range.==0]
     #c1_aff_equal = realfun_states[i][p_range.>0].== p_range[p_range.>0] 
     #c2_aff_equal = realfun_states[j][p_range.>0].== p_range[p_range.>0]

     c1_equal = realfun_states[i].== p_range 
     c2_equal = realfun_states[j].== p_range


     c1_higher_than_parent = all(!,realfun_states[i] .> p_range)
     c2_higher_than_parent = all(!,realfun_states[j] .> p_range)

     #check all zeros in parent range are zeros in daughter range (no jump dispersal)
     if all([c1_higher_than_parent, c2_higher_than_parent])

      #make sure all parental affinities are represented amongst daughter lineages (no affinites are lost totally via cladogeneis, just unequally inherited)
      if all((c2_equal .+ c1_equal) .> 0)
        #each daughter lineages has at least one real affinity from parent lineage 

        if in(2, p_range[c1_equal]) && in(2, p_range[c2_equal]) 





          push!(c1,realfun_states[i])
          push!(c2,realfun_states[j])
          print(i)
          print(" ")
          print(j)
          print("\n")


        end

      end

     end


     
    end  

  end 


  clados_unique=unique(eachrow((hcat([p_range for i in c1],c1,c2))))

  return(clados_unique)
end 







function get_clados(p_range::Array{Int64, 1}, realfun_states::Array{Array{Int64, 1}, 1},  retain_fun::Bool=true )                 

  #checks all biomes that the parent range is not realized in are the same, and the states within the parent realized rn ge are not zero

  if retain_fun
  index=[((p_range[p_range.!=2] == realfun_states[i][p_range.!=2]) && !in(0, realfun_states[i][p_range.==2]))  for i in 1:length(realfun_states)]

  else
  index=[((p_range[p_range.!=2] == realfun_states[i][p_range.!=2]) && !in(1, realfun_states[i][p_range.==2]))  for i in 1:length(realfun_states)]

  end

  nbiomes= length(realfun_states[1])


  if sum(p_range.==2)==1

    clados= hcat([p_range for i in 1],[p_range for i in 1], [p_range for i in 1])
    clados_filt=unique(eachrow(clados))

  else

    #generate all permutations of one child for a clado event

    c1=realfun_states[index]

    c1=[c1[i]  for i in 1:length(c1) 
      if c1[i]!=p_range]

    #realfun_states[index]

    #generate all corresponding children to c1
    c2=get_other_child_range(p_range, c1, retain_fun)


    #  get all possible cladogenetic splits structured as so vcat(hcat(parent range, all permutations of possible states, parent range)
    #                                                             hcat(parent range, parent range, all permutations of possible states)
    #                                                             hcat(parent range, all possible states, corresponding reverse states))
    #can add an additional column to specify par type (split verus partial full inheritance and all ecological versus vicariant )
    clados=vcat( hcat([p_range for i in c1],[p_range for i in c1], c1 ), 
                 hcat([p_range for i in c1], c1, [p_range for i in c1]), 
                 hcat([p_range for i in c1], c1 ,c2),
                 hcat([p_range for i in c1],[p_range for i in c1], [p_range for i in c1] )) 

    #filter to keep unique             
    clados_unique=unique(eachrow(clados))
    # remove all full sympatry for multiple realized biomes
    #clados_filt=[clados_unique[i]  for i in 1:length(clados_unique) 
    #  if !(clados_unique[i][1]==clados_unique[i][2] == clados_unique[i][3])]

        clados_filt=clados_unique
  end


  # add probabilities (currently uniform )
  #clados_filt=hcat(clados_filt,fill((1/length(realfun_states))/length(clados_filt),length(clados_filt) )) 


  #clados_filt=hcat(clados_filt,fill((1)/length(clados_filt),length(clados_filt) )) 


  return(clados_filt)
  #return(clados_unique)
end




function pstate_clados(p_range::Array{Int64, 1}, 
                     realfun_states::Array{Array{Int64, 1}, 1},
                     ecological,
                     allopatric)


  #filter to keep unique   
  fun_retain=[]

  if ecological
      push!(fun_retain,false)

  end
  if allopatric
      push!(fun_retain,true)

  end

  if !ecological && !allopatric

      return(print("error specify either ecological or allopatric cladogenesis"))

  end   

  print(p_range)
  clados=[]
  for clado_set in fun_retain

      #append!(clados, get_clados(p_range, realfun_states, clado_set))
      append!(clados, get_all_clados(p_range, realfun_states))


  end 

  clados_unique=unique(eachrow(clados))

  # add probabilities (currently uniform )
  #clados_filt=hcat(clados_filt,fill((1/length(realfun_states))/length(clados_filt),length(clados_filt) )) 


  clados_filt=hcat(clados_unique,fill((1)/length(clados_unique),length(clados_unique) )) 


  return(clados_filt)
end



"""
  make_state2slados(clado_Pmat::Matrix{Any})

TBW
"""
function make_state2clados(clado_Pmat::Matrix{Any})

  clados_by_P= Dict()

  for P_range in rf_states  

    clados_by_P[P_range]= [[clado_Pmat[i,1][1:end], clado_Pmat[i,2]]  for i in eachindex(clado_Pmat[:,1]) if clado_Pmat[i,1][1]==vec(P_range)]

  end

  return(clados_by_P)
end




function makeclado_Pmat_equal_prob(rf_states::Vector{Vector{Int64}},
                      ecological,
                      allopatric)

  rf2stDict = Dict(copy(rf_states).=>1:length(rf_states))

  cladoPmat=fill(0.0, length(rf_states),length(rf_states),length(rf_states))

  valid_clados=[pstate_clados(i,rf_states,ecological,allopatric ) for i in rf_states]


  for i in eachindex(valid_clados)
      jmax= Int(length(valid_clados[i])/2)
      for j in 1:jmax
          cladoPmat[rf2stDict[valid_clados[i][j][1][1]],rf2stDict[valid_clados[i][j][1][2]],rf2stDict[valid_clados[i][j][1][3]]]=valid_clados[i][(jmax)+j][1] 
      end
  end   
  return(cladoPmat)
end



"""
  get_clados_PMat(p_range::Array{Int64, 1}, 
                  realfun_states::Array{Array{Int64, 1}, 1})

Get full cladogenetic probability matrix for all possible cladogenetic events
given a specific parent range and generating all two child ranges 
"""
#function pstate_clados(p_range::Array{Int64, 1}, 
#                     realfun_states::Array{Array{Int64, 1}, 1})
#
#
#               
#
##checks all biomes that the parent range is n0t realized in are the same, and the states within the parent realized ra ge are not zero
#index=[((p_range[p_range.!=2] == realfun_states[i][p_range.!=2]) && !in(0, realfun_states[i][p_range.==2]))  for i in 1:length(realfun_states)]
#
#nbiomes= length(realfun_states[1])
#
#
#c1=realfun_states[index]
#
#if sum(p_range.==2)==1
#
#  clados= hcat([p_range for i in 1],[p_range for i in 1], [p_range for i in 1])
#  clados_filt=unique(eachrow(clados))
#
#else
#
#  #generate all permutations of one child for a clado event
#  c1=[c1[i]  for i in 1:length(c1) 
#    if c1[i]!=p_range]
#
#  #realfun_states[index]
#
#  #generate all corresponding children to c1
#  c2 =deepcopy(c1)
#
#  for i in eachindex(c1) for j in eachindex(c1[i])
#
#    #print(i)
#    #print(j)
#  
#      if c1[i][j]==2
#      
#        c2[i][j]=1
#      
#      elseif  c1[i][j]==1 && p_range[j]!=1
#      
#        c2[i][j]=2
#      
#       end 
#     
#  end  end
#
#  #  get all possible cladogenetic splits structured as so vcat(hcat(parent range, all permutations of possible states, parent range)
#  #                                                             hcat(parent range, parent range, all permutations of possible states)
#  #                                                             hcat(parent range, all possible states, corresponding reverse states))
#  clados=vcat( hcat([p_range for i in c1],[p_range for i in c1], c1 ), 
#               hcat([p_range for i in c1], c1, [p_range for i in c1]), 
#               hcat([p_range for i in c1], c1 ,c2),
#               hcat([p_range for i in c1],[p_range for i in c1], [p_range for i in c1] )) 
#
#  #filter to keep unique             
#  clados_unique=unique(eachrow(clados))
#  ## remove all full sympatry for multiple realized biomes
#  #clados_filt=[clados_unique[i]  for i in 1:length(clados_unique) 
#  #  if !(clados_unique[i][1]==clados_unique[i][2] == clados_unique[i][3])]
#
#
#end
#
#
## add probabilities (currently uniform )
##clados_filt=hcat(clados_filt,fill((1/length(realfun_states))/length(clados_filt),length(clados_filt) )) 
#
#
#clados_filt=hcat(clados_filt,fill((1)/length(clados_filt),length(clados_filt) )) 
#
#
#return(clados_filt)
#end




"""
  gen_cladoPMat(rf_states::Vector{Vector{Int64}})

TBW
"""
get_clados_sets(rf_states::Vector{Vector{Int64}}) =reduce( vcat, [pstate_clados(i,rf_states) for i in rf_states])



function fill_subsplit_cladoPmat(clado_Pmat, 
                               range_split_prob,
                               split_vec,
                               sub_vec)

  range_sub_prob=1-range_split_prob

  for p in eachindex(clado_Pmat[1,1,:])
      clado_slice=clado_Pmat[p,:,:]


      
      if length(sub_vec[p])>0 && length(split_vec[p])==0

        clado_slice[sub_vec[p]].= (1.0/length(sub_vec[p]))

        #clado_slice[split_vec[p]].=range_split_prob/length(split_vec[p])



      elseif length(sub_vec[p])==0 && length(split_vec[p])>0

        #clado_slice[sub_vec[p]].= (range_sub_prob/length(sub_vec[p]))

        clado_slice[split_vec[p]].= 1.0/length(split_vec[p])


      elseif length(sub_vec[p])>0 && length(split_vec[p])>0 

        clado_slice[sub_vec[p]]   .= range_sub_prob/length(sub_vec[p])

        clado_slice[split_vec[p]] .= range_split_prob/length(split_vec[p])

      end

      clado_Pmat[p,:,:]=clado_slice
  end

  return clado_Pmat

end

function fill_subsplit_cladoPmat(clado_Pmat, 
                                clado_probs,
                                split_vec,
                                sub_vec)



  range_split_prob = clado_probs[1]
  range_sub_prob   = clado_probs[2]
  range_equal_prob = clado_probs[3]

  for p in eachindex(clado_Pmat[1,1,:])
      clado_slice=clado_Pmat[p,:,:]


      
      if length(sub_vec[p])>0 && length(split_vec[p])==0

        #clado_slice[sub_vec[p]].= (1.0/length(sub_vec[p]))
        clado_slice[sub_vec[p]].= (range_sub_prob/length(sub_vec[p]))


        clado_slice[p,p]=range_equal_prob



      elseif length(sub_vec[p])==0 && length(split_vec[p])>0


        #clado_slice[split_vec[p]].= 1.0/length(split_vec[p])
        clado_slice[split_vec[p]].=range_split_prob/length(split_vec[p])

        clado_slice[p,p]=range_equal_prob



      elseif length(sub_vec[p])>0 && length(split_vec[p])>0 

        #print(sub_vec[p])
        #print((range_sub_prob))
        #print("\n")
        #print(length(sub_vec[p]))

        clado_slice[sub_vec[p]].= (range_sub_prob/length(sub_vec[p]))

        clado_slice[split_vec[p]].=range_split_prob/length(split_vec[p])

        clado_slice[p,p]=range_equal_prob


      end

      clado_Pmat[p,:,:]=clado_slice./sum(clado_slice)
  end

  return clado_Pmat

end


function make_clado_Pmat(rf_states, ecological, allopatric )


  clado_prior_vec=[Uniform(0,1),
                   Uniform(0,1)] 

  clado_prior_dists=make_Prior(clado_prior_vec, "split_sub_equal")

  #clado_probs_sim=[rand(clado_prior_dists[i]) for i in eachindex(clado_prior_dists)]
  clado_probs_sim=sample_sub_split_equal_clado_rates(clado_prior_dists)

  cladoPmat_unpar=makeclado_Pmat_equal_prob(rf_states, ecological, allopatric)

  String_Clado_Mats ,split_index_vec, sub_index_vec=get_cladoPmat_par_vecs(cladoPmat_unpar)

  print(clado_probs_sim[:,2])
  cladoPmat_sim=fill_subsplit_cladoPmat(cladoPmat_unpar, 
                                    clado_probs_sim[:,2],
                                    split_index_vec,
                                    sub_index_vec)


  return(cladoPmat_sim, clado_probs_sim, split_index_vec, sub_index_vec, String_Clado_Mats, clado_prior_dists)

end    

function sample_sub_split_equal_clado_rates(clado_prior_dists)
  clado_probs_sim=[0.0,0.0,0.0]
  #make sure clado_probs=1
  while (sum(clado_probs_sim)>1.0 || sum(clado_probs_sim)==0.0 )
    #split
    clado_probs_sim[1]=round(rand(clado_prior_dists[1]), digits=4)
    #sub
    clado_probs_sim[2]=round(rand(clado_prior_dists[2]), digits=4)
   # print(sum(clado_probs_sim))
  end

  
  clado_probs_sim[3]=round(1-clado_probs_sim[1]-clado_probs_sim[2], digits=4)

  clado_probs_sim=hcat(["split", "sub", "equal"], clado_probs_sim)

  return(clado_probs_sim)

end

function get_cladoPmat_par_vecs(clado_Pmat)
#creates arrays of cartesian indices for cladoPmat elements that are "cladogenetic range splits" or "cladogenetics range subset" 
split_vec  =  [CartesianIndex{2}[] for i in eachindex(clado_Pmat[:,1,1]) ]
sub_vec    =  [CartesianIndex{2}[] for i in eachindex(clado_Pmat[:,1,1]) ]

nstates=clado_Pmat[1,:,:][1,:]
String_Clado_Mats=fill("0", length(nstates),length(nstates),length(nstates))


for p in eachindex(clado_Pmat[:,1,1]) for lc in eachindex(clado_Pmat[1,:,1]) for rc in eachindex(clado_Pmat[1,1,:])
           if   clado_Pmat[p, lc, rc]>0
                if (p==lc==rc)
                  String_Clado_Mats[p, lc, rc]  =  "Equal"

                  
                elseif in(p, [lc,rc])
                    String_Clado_Mats[p, lc, rc]  =  "Sub"
                    #Par_Clado_Mats[p, lc, rc]     =    1
                    #push!(sub_vec[p], CartesianIndex(p, lc, rc))
                else 
                    String_Clado_Mats[p, lc, rc]  =  "Split"
                    #Par_Clado_Mats[p, lc, rc]     =     2 
                    #push!(split_vec[p], CartesianIndex(p, lc, rc))
                end 
            end    
        end
    end
    split_vec[p]  =  findall(x->x=="Split", String_Clado_Mats[p,:,:])
    sub_vec[p]    =  findall(x->x=="Sub"  , String_Clado_Mats[p,:,:])
end

return String_Clado_Mats ,split_vec, sub_vec

end




#######################

#function makeclado_Pmat(rf_states::Vector{Vector{Int64}},
#                        rf2stDict::Dict{Vector{Int64}, Int64})
#
#    cladoPmat=fill(0.0, length(rf_states),length(rf_states),length(rf_states))
#
#    valid_clados=[pstate_clados(i,rf_states) for i in rf_states]
#
#    for i in eachindex(valid_clados)
#        jmax= Int(length(valid_clados[i])/2)
#        for j in 1:jmax
#            cladoPmat[rf2stDict[valid_clados[i][j][1]],rf2stDict[valid_clados[i][j][2]],rf2stDict[valid_clados[i][j][3]]]=valid_clados[i][(jmax)+j] 
#        end
#    end   
#    return(cladoPmat)
#end

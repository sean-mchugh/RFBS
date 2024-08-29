#end



function sim_rf_tips(rf_states::Vector{Vector{Int64}},
                     Q:: Matrix{Float64},
                     cladoPmat::Array{Float64, 3},
                     tree_file     ::String,
                     ntips::Int64,
                     rescale_tree::Bool=true
            )


 
            
  if tree_file =="NA"

    tree, bts= make_ape_tree(ntips, 1.0,0.0)

    #tree, bts= make_ape_tree(ntips)


  else
          
    tree, bts = read_tree(tree_file)  
          
  end
  # sort according to branching times
  ed=tree.ed

  #writedlm("ets.txt",ets)
  #writedlm("bts.txt", brs)


  edges = cat(ed, [2*ntips ntips + 1], dims = 1)
  triads=maketriads(edges)
  #tip_biomes=tip_areas

  node_path  = get_trav_path_4prun(triads )

  
  br = branching_times(tree)

  # sort according to branching times
 
  
  brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
   

  #rescales the brs but not the tree itself since its
  if rescale_tree==true 

    brs[:,3:end]=brs[:,3:end]/maximum(brs[:,3:end])
   
  end
                     
  sim_node_path=reverse(node_path)
  #cladoPmat=makeclado_Pmat(rf_states,rf2stDict)
  ntips=Int(brs[1,1]-1)

  clado_events= fill("NaN", (2*(ntips)-1))
  edge_states= fill(0, (2*(ntips)-1))
  post_clado_states= fill(0, (2*(ntips)-1))


  states=collect(eachindex(rf_states))

  root=sample(states, Weights(fill(1.0, length(Q[1,:]))))
  #root=sample(states, Weights(exp(Q*1000)[1,:]))
  edge_states[ntips+1]=root

  #root_probs=exp(Q*1000)[:,1]
  #edge_states[1]=root


  for i in eachindex(sim_node_path)

    trav_node=sim_node_path[i]


    Px  = trav_node[1]

    #print(Px)
    #print("    ")
    #child nodes
    Cx1 = trav_node[2]
    Cx2 = trav_node[3]
    #reset likelihood for next node  

    #Ll     = zeros(nstates)
    #LR     = zeros(nstates) 
    #Pclado = zeros(nstates) 
    #get times between 
    tl = brs[brs[:,1].==Px .&& brs[:,2].==Cx1,3][1]
    tr = brs[brs[:,1].==Px .&& brs[:,2].==Cx2,3][1]

    #Ptleft=exp(Q*tl)
    #Ptright=exp(Q*tr)
    Ptleft= exponential!(Q*tl, ExpMethodHigham2005())
    Ptright= exponential!(Q*tr, ExpMethodHigham2005())




    possible_clados=findall(>(0), cladoPmat[edge_states[Px],:,:])


    child_states=sample(possible_clados, Weights(cladoPmat[edge_states[Px],:,:][possible_clados]))


    post_clado_states[trav_node[2]]=child_states[1]
    post_clado_states[trav_node[3]]=child_states[2]

    if in(edge_states[Px], [child_states[1], child_states[2]])

      clado_events[Px]="sub"

    else

      clado_events[Px]="split"

    end 


    edge_states[trav_node[2]] = sample(states, Weights(Ptleft[child_states[1],:]))


    edge_states[trav_node[3]] = sample(states, Weights( Ptright[child_states[2],:]))


    #cladoPmat[edge_states[Px],:,:]


  end

  return(edge_states[1:ntips], vcat(edge_states, post_clado_states), clado_events, tree, brs, node_path)

end

function sim_rf_tips_ana(rf_states::Vector{Vector{Int64}},
                     Q:: Matrix{Float64},
                     tree_file     ::String,
                     ntips::Int64,
                     rescale_tree::Bool=true
            )


 
            
  if tree_file =="NA"

    tree, bts= make_ape_tree(ntips, 1.0,0.0)
              
  else
          
    tree, bts = read_tree(tree_file)  
      
  end
  # sort according to branching times
  ed=tree.ed

  #writedlm("ets.txt",ets)
  #writedlm("bts.txt", brs)


  edges = cat(ed, [2*ntips ntips + 1], dims = 1)
  triads=maketriads(edges)
  #tip_biomes=tip_areas

  node_path  = get_trav_path_4prun(triads )

  
  br = branching_times(tree)
  # sort according to branching times
 
  
  brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
   

  #rescales the brs but not the tree itself since its
  if rescale_tree==true 

    brs[:,3:end]=brs[:,3:end]/maximum(brs[:,3:end])
   
  end
                     
  sim_node_path=reverse(node_path)
  #cladoPmat=makeclado_Pmat(rf_states,rf2stDict)
  ntips=Int(brs[1,1]-1)


  edge_states= fill(0, (2*(ntips)-1))


  states=collect(eachindex(rf_states))

  root=sample(states, Weights(exp(Q*10000)[1,:]))

  edge_states[ntips+1]=root

  #root_probs=exp(Q*1000)[1,:]
  #edge_states[1]=root


  for i in eachindex(sim_node_path)

    trav_node=sim_node_path[i]


    Px  = trav_node[1]

    #print(Px)
    #print("    ")
    #child nodes
    Cx1 = trav_node[2]
    Cx2 = trav_node[3]
    #reset likelihood for next node  

    #Ll     = zeros(nstates)
    #LR     = zeros(nstates) 
    #Pclado = zeros(nstates) 
    #get times between 
    tl = brs[brs[:,1].==Px .&& brs[:,2].==Cx1,3][1]
    tr = brs[brs[:,1].==Px .&& brs[:,2].==Cx2,3][1]

    Ptleft=exp(Q*tl)
    Ptright=exp(Q*tr)


    #possible_clados=findall(>(0), cladoPmat[edge_states[Px],:,:])

    
    #child_states=sample(possible_clados, Weights(cladoPmat[edge_states[Px],:,:][possible_clados]))



    edge_states[trav_node[2]] = sample(states, Weights(Ptleft[edge_states[Px],:]))


    edge_states[trav_node[3]] = sample(states, Weights( Ptright[edge_states[Px],:]))


    #cladoPmat[edge_states[Px],:,:]


  end

  return(edge_states[1:ntips], tree, brs, node_path)

end

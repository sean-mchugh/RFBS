

###pulled from TRIBE

#=

Tree utilities for ESSE

Ignacio Quintero Mächler

t(-_-t)

September 19 2017

=#



"""
Immutable type of an R tree `phylo` object type.
"""
struct rtree
  ed  ::Array{Int64,2}
  el  ::Array{Float64,1}
  tlab::Array{String,1}
  nnod::Int64
end





"""
    read_tree(tree_file::String; 
              order::String = "cladewise", 
              branching_times::Bool = true)

Function to read a tree using `RCall`
to call **ape** tree reading capabilities. 
"""
function read_tree(tree_file      ::String; 
                   order          ::String = "cladewise", 
                   branching_times::Bool = true)

  str = reval("""
                library(\"ape\")
                print('loaded ape')
                tree     <- read.nexus('$tree_file') 
                print('read tree')
                #tree     <- reorder(tree, order = '$order')
                print('ordered tree')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end






"""
    make_ape_tree(n::Int64, 
                  λ::Float64, 
                  μ::Float64; 
                  order::String = "cladewise", 
                  branching_times::Bool = true)

Make a phylogenetic tree using `phytools` in R. 
"""
function make_ape_tree(n              ::Int64, 
                       λ              ::Float64, 
                       μ              ::Float64; 
                       order          ::String = "cladewise", 
                       branching_times::Bool   = true)

  str = reval("""
                library(ape)
                library(phytools)
                tree     <- pbtree(n = $n, b = $λ, d = $μ, extant.only = TRUE)
                tree     <- reorder(tree, order = '$order')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end





"""
    maketriads(ed::Array{Int64,2})

Make edge triads given the tree. The first number is the parent, 
the second and third the children.
"""
function maketriads(ed::Array{Int64,2}; rev::Bool = false)

  # internal nodes
  ins = unique(ed[:,1])[1:(end-1)]::Array{Int64,1}

  rev && sort!(ins, rev = true)

  ed1 = ed[:,1]
  ed2 = ed[:,2]

  trios = Array{Int64,1}[]

  # for all internal nodes
  for i in ins
    daus = findall(isequal(i), ed1)
    pushfirst!(daus, findfirst(isequal(i), ed2))
    push!(trios, daus)
  end

  return trios::Array{Array{Int64,1},1}
end





"""
    abs_time_branches(el  ::Array{Float64,1}, 
                      ed  ::Array{Int64,2},
                      ntip::Int64)

Make array with absolute initial and end time for each 
branch. Time goes backwards with the present being `0.0`.
"""
function abs_time_branches(el  ::Array{Float64,1}, 
                           ed  ::Array{Int64,2},
                           ntip::Int64)

  @inbounds begin
    # make real time edges
    elrt = zeros(size(ed))

    # edges 1
    ed1 = ed[:,1]

    for i in axes(elrt,1)
      # is not tip
      if ed[i,2] > ntip
        elrt[i,2] = elrt[findfirst(isequal(ed[i,2]), ed1)::Int64,1]
        elrt[i,1] = elrt[i,2] + el[i]
      else
        elrt[i,1] = el[i]
      end
    end

  end

  return elrt
end





"""
    brts(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Get branching times for a tree in "cladewise" order. 
Time goes backwards with the present being `0.0`.
"""
function brts(el  ::Array{Float64,1}, 
              ed  ::Array{Int64,2},
              ntip::Int64)

  @inbounds begin
    # make real time edges

    # edges 1
    ed1 = ed[:,1]
    ed2 = ed[:,2]

    ne   = lastindex(ed1)
    nn   = zeros(ntip-1)
    intn = findall(map(x -> x > ntip, ed2))

    for i in intn
      nn[ed2[i]-ntip] = nn[ed1[i] - ntip] + el[i]
    end

    # tree height
    trh = nn[ed1[ne] - ntip] + el[ne]
    for i in Base.OneTo(ntip-1)
      nn[i] = trh - nn[i]
    end

  end

  return nn
end





"""
    tree_height(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Estimate tree height.
"""
function tree_height(el  ::Array{Float64,1}, 
                     ed  ::Array{Int64,2},
                     ntip::Int64)

  @inbounds begin
    # tree hight
    th  = 0.0
    ed2 = ed[:,2]

    da::Int64 = findfirst(isequal(1), ed2)

    # if the first branch reaches the tree height
    ed[da,1] == (ntip+1) && return el[da]

    while true
      th += el[da]

      pr = findfirst(isequal(ed[da,1]), ed2)::Int64

      if ed[pr,1] == (ntip+1)
        th += el[pr]
        break
      end

      da = findfirst(isequal(ed[pr,2]), ed2)::Int64
    end
  end

  return th
end





"""
    postorderedges(ed  ::Array{Int64,2},
                   el  ::Array{Float64,1},
                   ntip::Int64)

Organize edges, edge lengths in postorder traversal fashion.
"""
function postorderedges(ed  ::Array{Int64,2},
                        el  ::Array{Float64,1},
                        ntip::Int64)

  # post-order transversal using 2 stacks
  stack1 = [ed[1]]
  stack2 = Int64[]

  while lastindex(stack1) > 0
    nod = pop!(stack1)
    push!(stack2, nod)

    if nod <= ntip
      continue
    else
      wn = findfirst(isequal(nod), ed[:,1])
      push!(stack1, ed[wn,2],ed[(wn+1),2])
    end
  end

  # rearrange edges accordingly
  indx = deleteat!(indexin(reverse(stack2), ed[:,2]),size(ed,1)+1)

  ed = ed[indx,:]
  el = el[indx]

  # advance nodes with only daughter tips
  tnd = Int64[]
  ndp = Int64[]
  for nd in unique(ed[:,1])
    fed = findall(ed[:,1 ] .== nd)
    if length(filter(x -> x <= ntip, ed[fed,2])) == 2
      push!(tnd, nd)
      push!(ndp, fed[1], fed[2])
    end
  end

  append!(ndp,setdiff(1:(2ntip-2), ndp))

  ed[ndp,:]

  return ed[ndp,:], el[ndp]
end







"""
    make_ape_tree(n::Int64, 
                  λ::Float64, 
                  μ::Float64; 
                  order::String = "cladewise", 
                  branching_times::Bool = true)

Make a phylogenetic tree using `ape` in R. 
"""
function make_ape_tree_RIS(n              ::Int64, 
                       order          ::String = "cladewise", 
                       branching_times::Bool   = true)

  str = reval("""
                library(ape)
                tree     <- ape::chronos(rtree(n='$n'))
                tree     <- reorder(tree, order = '$order')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end



"""
    make_ape_tree(n::Int64, 
                  λ::Float64, 
                  μ::Float64; 
                  order::String = "cladewise", 
                  branching_times::Bool = true)

Make a phylogenetic tree using `phytools` in R. (use when not on RIS phytools doesnt work on R version 4.0.4 easily)
"""
function make_ape_tree(n              ::Int64, 
                       λ              ::Float64, 
                       μ              ::Float64; 
                       order          ::String = "cladewise", 
                       branching_times::Bool   = true)

  str = reval("""
                library(ape)
              library(phytools)
               tree     <- pbtree(n = $n, b = $λ, d = $μ, extant.only = TRUE)
                tree     <- reorder(tree, order = '$order')
                edge     <- .subset2(tree,'edge')
                Nnode    <- .subset2(tree,'Nnode')
                tiplabel <- .subset2(tree,'tip.label')
                edlength <- .subset2(tree,'edge.length')
                list(edge,Nnode,tiplabel,edlength)
              """)

  edge     = rcopy(str[1])
  edge     = convert(Array{Int64},edge)
  Nnode    = rcopy(str[2])
  Nnode    = convert(Int64,Nnode)
  tiplabel = rcopy(str[3])
  edlength = rcopy(str[4])
  edlength = convert(Array{Float64},edlength)

  tree = rtree(edge, edlength, tiplabel, Nnode)

  if branching_times
    brtimes = reval("""
                      brtimes <- branching.times(tree)
                    """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
  else
    return tree
  end
end



"""
    maketriads(ed::Array{Int64,2})

Make edge triads given the tree. The first number is the parent, 
the second and third the children.
"""
function maketriads(ed::Array{Int64,2}; rev::Bool = false)

  # internal nodes
  ins = unique(ed[:,1])[1:(end-1)]::Array{Int64,1}

  rev && sort!(ins, rev = true)

  ed1 = ed[:,1]
  ed2 = ed[:,2]

  trios = Array{Int64,1}[]

  # for all internal nodes
  for i in ins
    daus = append!([i],ed2[ed1.==i])
    
    push!(trios, daus)
  end

  return trios::Array{Array{Int64,1},1}
end




"""
    get_trav_path_4prun(tip_biomes::Dict{Int64, Vector{Int64}},
                             rf2stDict ::Dict{Vector{Int64}, Int64},
                             triads    ::Vector{Vector{Int64}}
                             )

TBW
"""
function get_trav_path_4prun(
                             triads    ::Vector{Vector{Int64}}
                             )
  ntips=length(triads)+1
  #nstates=length(rf2stDict)

      #setup initial tip state probabilities multiplying by 1.0 to make floats
  #edge_end_prob=Dict{Int64, Any}((1:length(tip_biomes)) .=>  [ 1.0*((1:length(rf2stDict)).==rf2stDict[tip_biomes[i]]) for i in 1:length(tip_biomes)])

  covered_nodes=collect(1:ntips)

  covered_triads=[triads[j] for j in 1:length(triads) if (in(triads[j][2], covered_nodes) .&& in(triads[j][3], covered_nodes) .&& !in(triads[j][1], covered_nodes) )]


  

  while length(covered_nodes) < (ntips*2-1) 



  #identify set of traversable nodes, starting nodes connecting to only tips and then nodes 
  #connecting to the nodes with node state probs or tips
    trav_nodes=[triads[j] for j in 1:length(triads) if (in(triads[j][2], covered_nodes)
    .&& in(triads[j][3], covered_nodes)
    .&& !in(triads[j][1], covered_nodes) )]

    for node in eachindex(trav_nodes)

      push!(covered_nodes,trav_nodes[node][1])
      push!(covered_triads,trav_nodes[node])


    end

  end

  return(unique(covered_triads))
end




"""
    get_clados_PMat(p_range::Array{Int64, 1}, 
                    realfun_states::Array{Array{Int64, 1}, 1})

Get full cladogenetic probability matrix for all possible cladogenetic events
 given a specific parent range and generating all two child ranges 
"""
function pstate_clados(p_range::Array{Int64, 1}, 
                       realfun_states::Array{Array{Int64, 1}, 1})


                 

  index=[((p_range[p_range.!=2] == realfun_states[i][p_range.!=2]) && !in(0, realfun_states[i][p_range.==2]))  for i in 1:length(realfun_states)]

  nbiomes= length(realfun_states[1])


  c1=realfun_states[index]

  if sum(p_range.==2)==1

    clados= hcat([p_range for i in 1],[p_range for i in 1], [p_range for i in 1])
    clados_filt=unique(eachrow(clados))



  else

    c1=[c1[i]  for i in 1:length(c1) 
      if c1[i]!=p_range]

    #realfun_states[index]

    c2 =deepcopy(c1)

    for i in eachindex(c1) for j in eachindex(c1[i])

      #print(i)
      #print(j)
    
        if c1[i][j]==2
        
          c2[i][j]=1
        
        elseif  c1[i][j]==1 && p_range[j]!=1
        
          c2[i][j]=2
        
         end 
       
    end  end

    #  get all possible cladogenetic splits structured as so vcat(hcat(parent range, all permutations of possible states, parent range)
    #                                                             hcat(parent range, parent range, all permutations of possible states)
    #                                                             hcat(parent range, all possible states, corresponding reverse states))
    clados=vcat( hcat([p_range for i in c1],[p_range for i in c1], c1 ), 
                 hcat([p_range for i in c1], c1, [p_range for i in c1]), 
                 hcat([p_range for i in c1], c1 ,c2)) 

    #filter to keep unique             
    clados_unique=unique(eachrow(clados))
    # remove all full sympatry for multiple realized biomes
    clados_filt=[clados_unique[i]  for i in 1:length(clados_unique) 
      if !(clados_unique[i][1]==clados_unique[i][2] == clados_unique[i][3])]


  end


  # add probabilities (currently uniform )
  #clados_filt=hcat(clados_filt,fill((1/length(realfun_states))/length(clados_filt),length(clados_filt) )) 


  clados_filt=hcat(clados_filt,fill((1)/length(clados_filt),length(clados_filt) )) 


  return(clados_filt)
end


"""
    branching_times(tree::rtree)

Function to estimate absolute
branching times, with time 0 at the
present, time going backward.
"""
function branching_times(tree::rtree)

  @inbounds begin

    n    = tree.nnod + 1
    el_t = findall(x -> x <= n, tree.ed[:,2])

    brs = zeros(lastindex(tree.el),5)

    brs[:, 1:2] = tree.ed
    brs[:, 3]   = tree.el

    for j = eachindex(tree.el)
      if brs[j,1] == (n+1)
        brs[j,4] = 0.0
        brs[j,5] = brs[j,4] + brs[j,3]
      else
        brs[j,4] = brs[brs[j,1] .== brs[:,2],5][1]
        brs[j,5] = brs[j,4] + brs[j,3]
      end
    end

     # change time forward order
    @views tree_depth = brs[n .== brs[:,2],5][1]

    for j = eachindex(tree.el) 
      brs[j,4] = tree_depth - brs[j,4]
      brs[j,5] = tree_depth - brs[j,5]
    end

    brs[el_t,5] .= 0.0

  end

  return brs
end

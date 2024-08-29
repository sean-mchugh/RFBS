#function rfprune(rate_pars        ,
#    Q_par_matrix      ,
#    Q_mat_empty      ,
#    Q_par_index_vec  ,
#    tip_probs        ,
#    clado_Pmat       ,
#    brs              ,
#    node_path        
#   )
#
#
#
##print(Q_rates)
#
##map rates to Q matrix 
#Q= fill_emptyQmat(Q_mat_empty, rate_pars, Q_par_index_vec, Q_par_matrix)
#
##print(Q)
#
#lL, node_probs, post_clado_probs = rf_prune_ll_new(Q ,          
#          tip_probs,  
#          clado_Pmat,  
#          brs,         
#          node_path )
#
#
#
#
#return(lL, node_probs,post_clado_probs)
#
#
#end


function open_new_log(log_filename)
    if isfile(log_filename*".txt")
        i=1

        while isfile(log_filename*string(i)*".txt")

            i+=1

        end

        log_filename_txt=log_filename*string(i)*".txt"
        open(log_filename_txt, "w") do nadda

        end   
    else

        log_filename_txt=log_filename*".txt"
        open(log_filename_txt, "w") do nadda

        end   


    end
    return(log_filename_txt)

end


#cladogenetic likelihood calvulcation for single cladogenesis event on phylogeny
function Cpartial_ll(left_child_states::Array{Float64, 1},
                     right_child_states::Array{Float64, 1},
                          cladoPmat::Array{Float64, 3})


    parent_states=fill(0.0, length(left_child_states))

    for i in eachindex(parent_states) for j in eachindex(left_child_states) for z in eachindex(right_child_states)

                parent_states[i]+= cladoPmat[i,j,z]*left_child_states[j]*right_child_states[z]

            end

        end
    end
    return( parent_states::Vector{Float64})
end


#anagenetic likelihood calculation for single branch on phylogeny
function Apartial_ll(end_states::Array{Float64, 1},
                     QPmat::Array{Float64,2})


    start_states=fill(0.0, length(end_states))

    for i in eachindex(start_states) for j in eachindex(end_states)

            start_states[i]+= QPmat[i,j]*end_states[j]

        end

    end
    return(start_states::Vector{Float64})
end


function make_Prior(prior_vec, move_types)
    if length(prior_vec)>1

        ratePrior=[ par for par in prior_vec]

    else

        ratePrior=[ prior_vec[1] for par in eachindex(move_types)]

    end

    return(ratePrior)

end



#Q=Q_start      
#tip_probs    
#cladoPmat =clado_Pmat_c 
#brs          
#node_path    

function rf_prune_algo(Q          ::  Matrix{Float64},
                    tip_probs       ::  Vector{Vector{Float64}},
                    cladoPmat       ::  Array{Float64, 3},
                    #triads         ::  Vector{Vector{Int64}},
                    brs            ::  Matrix{Float64},
                    node_path      ::Vector{Vector{Int64}}
                    #par2clado      ::Dict{Any, Any}
                    )
  
    #edge_end_prob=Dict{Int64, Any}((1:length(tip_biomes)) .=>  [ 1.0*((1:length(rf2stDict)).==rf2stDict[tip_biomes[i]]) for i in 1:length(tip_biomes)])


    uf_thresh=10^-100

    uf_counter=0

    ntips=length(tip_probs)

    nstates=length(Q[1,:])

    edge_probs=vcat(tip_probs, [fill(0.0, length(Q[1,:])) for i in 1:(length(tip_probs)-1) ])
    post_clado_probs=vcat(tip_probs, [fill(0.0, length(Q[1,:])) for i in 1:(length(tip_probs)-1) ])

    #get stationary probs by getting the Prob matrix over a long period of time
    #root_probs=exp(Q*1000)[1,:]
    root_probs=fill(1.0, length(Q[1,:]))



    logL=zeros(1)

    #cladoPmat=makeclado_Pmat(rf_states,rf2stDict)

    for i in eachindex(node_path)

        trav_node=node_path[i]

        Px  = trav_node[1]

        #print(Px)
        #print("    ")
        #child nodes
        Cx1 = trav_node[2]
        Cx2 = trav_node[3]
        #reset likelihood for next node  

        Ll     = zeros(nstates)
        LR     = zeros(nstates) 
        Pclado = zeros(nstates) 
        #get times between 
        tl = brs[brs[:,1].==Px .&& brs[:,2].==Cx1,3][1]
        tr = brs[brs[:,1].==Px .&& brs[:,2].==Cx2,3][1]

        #Pt_l=exp(Q*tl)
        #Pt_r=exp(Q*tr)

        Pt_l= exponential!(Q*tl, ExpMethodHigham2005())
        Pt_r= exponential!(Q*tr, ExpMethodHigham2005())

       
        edge_probs[Cx1].<uf_thresh

        edge_probs[Cx2]

        Ll=Apartial_ll(edge_probs[Cx1], Pt_l)

        if sum(Ll.<uf_thresh)>0
            Ll=Ll.*(10.0^100.0)
            uf_counter= 1+uf_counter
        end    



        Lr=Apartial_ll(edge_probs[Cx2], Pt_r)

        if sum(Lr.<uf_thresh)>0
            Lr=Lr.*(10.0^100.0)
            uf_counter= 1+uf_counter
        end    

        post_clado_probs[Cx2]=Lr
        post_clado_probs[Cx1]=Ll
    

        edge_probs[Px]=Cpartial_ll(Ll,Lr,cladoPmat)

        if(sum((edge_probs[Px]).==Inf)>0) 
            #print(i)
            #print("\n")
            #print(edge_probs[Px])
        end


        if(sum(isnan.(edge_probs[Px]))>0) 
            #print(i)
            #print("\n")
           # print(edge_probs[Px])
        end
    end


    #print(edge_probs[end])

    logL=log(sum(root_probs.*edge_probs[ntips+1]))+(-100*uf_counter*log(10))
    #logL=sum(log.(root_probs.*edge_probs[ntips+1]))+(-100*uf_counter*log(10))
    #logL=sum(log.(sum.(edge_probs)))+(-100*uf_counter*log(10))

    return(logL::Float64, edge_probs::Vector{Vector{Float64}},post_clado_probs::Vector{Vector{Float64}} )

end



function sample_anc_states(edge_probs,  post_clado_probs, root_probs, tip_probs, Q_mat_empty, rate_pars, Q_par_index_vec, Q_par_matrix, cladoPmat, brs, node_path)


    Q= fill_emptyQmat(Q_mat_empty, rate_pars, Q_par_index_vec, Q_par_matrix)


    nstates=length(Q[1,:])
    ntips=Integer((length(edge_probs)+1)/2)
    root_final_probs=root_probs .* edge_probs[ntips+1]

    #print(edge_probs[ntips+1])
    #print("\n") 
 #
    #print(root_final_probs)
    #print("\n") 


    anc_root_probs_sum=sum(root_final_probs)
    anc_root_probs= fill(0.0, size(root_probs))


    clado_events= fill("NaN", length(edge_probs))

    Nst= fill(0, length(edge_probs) )
    Cst= fill(0, length(edge_probs) )
    
    #print(rf_states)

    for i in eachindex(root_final_probs) 



      anc_root_probs[i]= root_final_probs[i]/anc_root_probs_sum

    end 



    #print(anc_root_probs)
    #print("\n") 

    #sample node state for root

    Nst[ntips+1]=sample(eachindex(anc_root_probs),Weights(anc_root_probs))

    #generate traversale path fro ancestral state recon
    anc_st_node_path=reverse(node_path)

    #generate all possible clados as matrix of children state states and probs with each row being a clado scenario
    clado_probs_vecs=[hcat(findall(item -> item>0, cladoPmat[state,:,:]), cladoPmat[state,findall(item -> item>0, cladoPmat[state,:,:])]) for state in eachindex(root_final_probs)]


    for i in eachindex(anc_st_node_path)

        trav_node = anc_st_node_path[i]

        Px  = trav_node[1]
        Pst = Nst[Px]
        #print(Px)
        #print("    ")
        #child nodes
        Cx1 = trav_node[2]
        Cx2 = trav_node[3]
        #reset likelihood for next node  

        Ll     = zeros(nstates)
        LR     = zeros(nstates) 
        #Pclado = zeros(nstates) 
        #get times between 
        tl = brs[brs[:,1].==Px .&& brs[:,2].==Cx1,3][1]
        tr = brs[brs[:,1].==Px .&& brs[:,2].==Cx2,3][1]

        Ptl=exp(Q*tl)
        Ptr=exp(Q*tr)


        #draw a cladogenetic scenario based on the drawn parent   

        anc_clado_probs=copy(clado_probs_vecs[Pst])
        for clado in 1:size(clado_probs_vecs[Pst])[1]


            clado_prob=clado_probs_vecs[Pst][clado,:]
            anc_clado_probs[clado,2]=clado_prob[2]*post_clado_probs[Cx1][clado_prob[1][1]]*post_clado_probs[Cx2][clado_prob[1][2]]

        end

        anc_clado_probs_sum=sum(anc_clado_probs[:,2])
        anc_clado_probs_bysum=copy(clado_probs_vecs[Pst])

        for clado in 1:size(clado_probs_vecs[Pst])[1]


            anc_clado_probs_bysum[clado,2]=anc_clado_probs[clado,2]/anc_clado_probs_sum

        end

        anc_clado=sample(clado_probs_vecs[Pst][:,1], Weights(Float64.(anc_clado_probs_bysum[:,2])))


        Cst[Cx1]=anc_clado[1]
        Cst[Cx2]=anc_clado[2]


        if in(Pst, [Cst[Cx1],Cst[Cx2]])

            if Pst==Cst[Cx1]==Cst[Cx2]

                clado_events[Px]="equal"

            else   
                clado_events[Px]="sub"
      
            end
          else
      
            clado_events[Px]="split"
      
          end 


        anc_ana_probsCx1=zeros(eachindex(root_final_probs))
        anc_ana_probsCx2=zeros(eachindex(root_final_probs))

        for ana in eachindex(root_final_probs)

            anc_ana_probsCx1[ana]= Ptl[ Cst[Cx1] ,ana]*edge_probs[Cx1][ana]
            anc_ana_probsCx2[ana]=Ptr[ Cst[Cx2], ana]*edge_probs[Cx2][ana]

        end

        anc_ana_probs_sumCx1=sum( anc_ana_probsCx1)
        anc_ana_probs_bysumCx1=zeros(eachindex(root_final_probs))

        anc_ana_probs_sumCx2=sum( anc_ana_probsCx2)
        anc_ana_probs_bysumCx2=zeros(eachindex(root_final_probs))


        for ana in eachindex(root_final_probs)


            anc_ana_probs_bysumCx1[ana]= anc_ana_probsCx1[ana]/anc_ana_probs_sumCx1
            anc_ana_probs_bysumCx2[ana]= anc_ana_probsCx2[ana]/anc_ana_probs_sumCx2

        end

        if in(Cx1, 1:ntips)


            Nst[Cx1]=sample( eachindex(anc_ana_probs_bysumCx1), Weights(Float64.(anc_ana_probs_bysumCx1) .* tip_probs[Cx1]))


        else 

            Nst[Cx1]=sample( eachindex(anc_ana_probs_bysumCx1), Weights(Float64.(anc_ana_probs_bysumCx1)  ))
        end



        if in(Cx2, 1:ntips)


            Nst[Cx2]=sample( eachindex(anc_ana_probs_bysumCx2), Weights(Float64.(anc_ana_probs_bysumCx2) .* tip_probs[Cx2]))

        else 

            Nst[Cx2]=sample( eachindex(anc_ana_probs_bysumCx2), Weights(Float64.(anc_ana_probs_bysumCx2)))
        end


      


    end


    return(vcat(Nst,Cst), clado_events)


end 






function   rate_par_upd(rate_pars_c      ,
                        tip_probs        ,
                        clado_Pmat_c     ,
                        lL_c             ,
                        prunelL_c        ,
                        clado_priorlL_c  ,
                        rate_priorlL_c   , 
                        edge_probs  ,
                        post_clado_probs,                      
                        Q_par_matrix     ,
                        Q_zeros          ,
                        Q_index_vec  ,
                        prior_dists      ,
                        brs              ,
                        node_path        ,
                        tuning_par_vec   ,
                        par_acceptfreq_vec   ,
                        par_propfreq_vec  ,
                        prior_only=false
                        )

    rate_pars_p=copy(rate_pars_c)
    lL_p       =copy(lL_c       )
    prunelL_p  =copy(prunelL_c  )
    rate_priorlL_p  =copy(rate_priorlL_c  )


    prop_par_ind=rand(eachindex(rate_pars_c))

    rate_pars_p[prop_par_ind], hr   = multi_move(rate_pars_c[prop_par_ind], #parameter
                                                 tuning_par_vec[prop_par_ind]  #tuning par
                                            )
    


     Q_p= fill_emptyQmat( Q_zeros, rate_pars_p, Q_index_vec, Q_par_matrix)

     prunelL_p, edge_probs_p, post_clado_probs_p = rf_prune_algo(Q_p     ,
                                                     tip_probs        ,
                                                     clado_Pmat_c       ,
                                                     brs              ,
                                                     node_path        
                                                              )
    
    
    rate_priorlL_p = sum([log(pdf(prior_dists[i], rate_pars_p[i])) for i in eachindex(rate_pars_p)]) 


    if prior_only 

        lL_p= rate_priorlL_p + clado_priorlL_c +hr

    else

        lL_p= rate_priorlL_p + clado_priorlL_c + prunelL_p  +hr

    end 


    lL_ratio =  exp(lL_p-lL_c)

    if rand(Uniform(0, 1)) < lL_ratio 

        #Q_c         =copy(Q_p)
        rate_pars_c  =copy(rate_pars_p)
        rate_priorlL_c    =copy(rate_priorlL_p )
        prunelL_c  =copy(prunelL_p  )
        lL_c         =copy(lL_p )
        edge_probs      = edge_probs_p         
        post_clado_probs= post_clado_probs_p 




        par_acceptfreq_vec[prop_par_ind] = par_acceptfreq_vec[prop_par_ind] .+1
    end 

    par_propfreq_vec[prop_par_ind] = par_propfreq_vec[prop_par_ind] .+1

    return (
            rate_pars_c ,
            edge_probs  ,
            post_clado_probs,
            rate_priorlL_c   ,
            prunelL_c   , 
            lL_c        ,
            par_acceptfreq_vec,
            par_propfreq_vec  )

end

#clado_Pmat_unpar=cladoPmat_unpar




function   clado_par_upd(rate_pars_c      ,
                         clado_probs_c     , 
                         clado_Pmat_c           ,
                         lL_c                   ,
                         prunelL_c              ,
                         clado_priorlL_c        ,
                         rate_priorlL_c         ,  
                         edge_probs         ,
                         post_clado_probs   ,
                         clado_prior            , 
                         split_vec              ,
                         sub_vec                ,
                         Q_par_matrix     ,
                         Q_zeros          ,
                         Q_index_vec  ,
                         tip_probs              ,
                         brs                    ,
                         node_path              ,
                         clado_tuning_par_vec   ,
                         par_acceptfreq_vec     ,
                         par_propfreq_vec       ,
                         prior_only=false
                         )

    clado_probs_p     =copy(clado_probs_c )
    lL_p             =copy(lL_c       )
    prunelL_p        =copy(prunelL_c  )
    clado_priorlL_p  =copy(clado_priorlL_c  )
    #rate_priorlL_p  =copy(rate_priorlL_c  )                      


    #prop_par_ind=rand(eachindex(rate_pars_c))

    #rate_pars_p[prop_par_ind], hr   =     multi_move( rate_pars_c[prop_par_ind], #parameter
    #                                                    tuning_par_vec[prop_par_ind]  #tuning par
    #                                                  )
    #
     #print(Q)

     #choose which clado event prob to propose on (equal will always update with either split or sub)

     par_up=rand(1:2)

     #print(clado_probs_c[par_up])
     clado_probs_p[par_up], hr     = clado_multi_move(clado_probs_c[par_up], clado_tuning_par_vec[1])
    
     clado_probs_p[3]=1-sum(clado_probs_p[1:2])

    if  0 < sum(clado_probs_p[1:2]) < 1.0 

        Q= fill_emptyQmat( Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)


         clado_Pmat_p            = fill_subsplit_cladoPmat(clado_Pmat_c, 
                                                           clado_probs_p,
                                                           split_vec,
                                                           sub_vec)

         prunelL_p, edge_probs_p, post_clado_probs_p = rf_prune_algo(Q     ,
                                                                tip_probs        ,
                                                                clado_Pmat_p       ,
                                                                brs              ,
                                                                node_path        
                                                                         )
        

        
        clado_priorlL_p = sum([log(pdf(clado_prior[i], clado_probs_p[i])) for i in eachindex(clado_prior)]) 


        lL_p= clado_priorlL_p + rate_priorlL_c + prunelL_p + hr

        if prior_only 

            lL_p=clado_priorlL_p + rate_priorlL_c + hr

        else

            lL_p=clado_priorlL_p + rate_priorlL_c + prunelL_p + hr

        end 


        lL_ratio =  exp(lL_p-lL_c)

        if rand(Uniform(0, 1)) < lL_ratio 

            clado_probs_c     =copy(clado_probs_p   )
            clado_priorlL_c  =copy(clado_priorlL_p)
            prunelL_c        =copy(prunelL_p      )
            lL_c             =copy(lL_p           )
            edge_probs      = edge_probs_p         
            post_clado_probs= post_clado_probs_p 

            par_acceptfreq_vec[[par_up,3]] = par_acceptfreq_vec[[par_up,3]] .+1
        
        
        end 

    end

    par_propfreq_vec[[par_up,3]] = par_propfreq_vec[[par_up,3]] .+1

    return (clado_probs_c       ,
            edge_probs         ,
            post_clado_probs   ,
            clado_priorlL_c          ,
            prunelL_c    , 
            lL_c               ,
            par_acceptfreq_vec ,
            par_propfreq_vec   )

end


#clado_Pmat_unpar= cladoPmat_unpar
#start_clado_split_prob=.2


function realfun_mcmc(start_rates, start_clado_probs, 
                      iters, iter_trims, write_interval,
                      anc_state_sampling, 
                      Q_par_matrix, Q_index_vec, move_types, rate_prior_dists, 
                      clado_Pmat_unpar,  clado_prior_dists,
                      rf_states,
                      tip_probs, 
                      tree, 
                      log_filename, 
                      rate_tuning_par, clado_tuning_par,
                      rate_par_proposal_prob,
                      prior_only=false,
                      rescale_tree=true)



    #make prior distribution vectors for rates a cladogenetic splitting probability 
    #rate_prior_dists=make_Prior(rate_prior_vec, move_types)
    #clado_prior_dists=make_Prior(clado_prior_vec, clado_types)
    #print(tree)

    #tip_probs=get_tip_probs(rf_states,tip_states)

    ed=tree.ed
    ntips=length(tree.tlab)
    #writedlm("ets.txt",ets)
    #writedlm("bts.txt", brs)

   # clado_Pmat=makeclado_Pmat(rf_states,rf2stDict)


    edges = cat(ed, [2*ntips ntips + 1], dims = 1)
    triads=maketriads(edges)
    #tip_biomes=tip_areas
    node_path  = get_trav_path_4prun(triads )
    br = branching_times(tree)
    # sort according to branching times
    brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)
    
    if(rescale_tree)
        brs[:,3:end]=brs[:,3:end]/maximum(brs[:,3:end])
    end
    root_probs= fill(1.0, length(rf_states))
    #lower = optim_lower_bound
    #upper = optim_upper_bound
    #initial_x = [0.2]
    #inner_optimizer = GradientDescent()
    #start_rates =[rand(prior_dists[i]) for i in eachindex(move_types)]
    #start_rates =[1.0 for i in eachindex(move_types)]



    #generate current Q matrix (doing it at mcmc level instead of ll fn level so it doesnt need to be reassembled every time a par is updated (say a cladogenetic par))
    Q_zeros =zeros(length(rf_states), length(rf_states))
    rate_pars_c=copy(start_rates)
    Q_start= fill_emptyQmat( Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix)

    ##generate cladogenetic probability matrix and set as current matrix (either (1) use an unparameterized matrix where all daughter scenarios are equal for any one given parental state, or (2) one where it is weighted by range splitting vs range subsetting

    #clado_split_prob_c=copy(start_clado_split_prob)
    #clado_sub_prob_c=copy(start_clado_sub_prob)
    clado_probs_c=copy(start_clado_probs[:,2])
    String_Clado_Mats, split_index_vec, sub_index_vec=get_cladoPmat_par_vecs(clado_Pmat_unpar)
    #print(clado_probs_c)
    clado_Pmat_c=fill_subsplit_cladoPmat(clado_Pmat_unpar, 
                                         start_clado_probs[:,2],
                                         split_index_vec,
                                         sub_index_vec)


    ###generate tuning parameterrs, aceptance frequencies and proposal frequencies for parameters
    rate_tuning_par_vec= [rate_tuning_par for i in move_types]
    rate_par_acceptfreq_vec= [0 for i in rate_pars_c]
    rate_par_propfreq_vec= [0 for i in rate_pars_c]


    clado_tuning_par_vec          =  [clado_tuning_par for i in clado_probs_c]
    clado_par_acceptfreq_vec      =  [0 for i in clado_probs_c]
    clado_par_propfreq_vec        =  [0 for i in clado_probs_c]

    AR_move_types=["AR"*move for move in move_types]
    #log_header_string= "iter" * "\t" * "lL_c" * "\t" * "prunelL_c" * "\t" * "priorlL_c" * "\t" *join( move_types, "\t") * "\t"  * join(  AR_move_types, "\t")
    log_header_string= "iter" * "\t" * 
                        "lL_c" * "\t" *
                        "prunelL_c" * "\t" * 
                        "rate_priorlL_c" * "\t" * 
                        "clado_priorlL_c" * "\t" * 
                        join( move_types, "\t") * "\t"  * 
                        join( start_clado_probs[:,1], "\t") * "\t"  * 
                        join(  AR_move_types, "\t") *"\t" *
                        "AR_clado_split_prob" *"\t" *
                        "AR_clado_sub_prob" *"\t" *
                        "AR_clado_equal_prob" 



    iter=0

    file_line=log_header_string

    anc_states_log_filename=log_filename*"_anc_states"
    anc_clado_log_filename=log_filename*"_anc_clados"

    #log_filename = open_new_log(log_filename)
    #anc_states_log_filename = open_new_log( anc_st_log_filename)
    #anc_clado_log_filename = open_new_log( anc_clado_log_filename)


    chain_cache=[]
    anc_clado_cache=[]
    anc_states_cache=[]

    #open(log_filename, "w") do log_file

    #    open(anc_states_log_filename, "w") do anc_state_log_file

    #        open(anc_clado_log_filename, "w") do anc_clado_log_file



        log_file =  open(log_filename,"a")
        #write(exampleFileIOStream, file_line)
        println(log_file, file_line)
        print(file_line)

        close(log_file)


    




        prunelL_c, edge_probs,  post_clado_probs = rf_prune_algo(Q_start      ,
                                                                 tip_probs    ,
                                                                 clado_Pmat_c ,
                                                                 brs          ,
                                                                 node_path        
                                                                         )


        rate_priorlL_c = sum([log(pdf(rate_prior_dists[i], rate_pars_c[i])) for i in eachindex(rate_pars_c)]) 

        clado_priorlL_c = sum([log(pdf(clado_prior_dists[i],clado_probs_c[i])) for i in eachindex(clado_prior_dists)]) 

        if prior_only 

            lL_c= rate_priorlL_c + clado_priorlL_c

        else

            lL_c= rate_priorlL_c + clado_priorlL_c + prunelL_c  

        end 
        log_line_string= "0" * "\t" * 
                         string(round(lL_c, digits=4)) * "\t" * 
                         string(round(prunelL_c,digits=4)) * "\t" * 
                         string(round(rate_priorlL_c,digits=4)) * "\t" * 
                         string(round(clado_priorlL_c, digits=4)) * "\t" *
                         join( string.(round.(rate_pars_c, digits=5)), "\t") * "\t"  * 
                         join( string.(round.(clado_probs_c, digits=5)), "\t") * "\t"  * 
                         join( string.(round.(rate_par_acceptfreq_vec ./rate_par_propfreq_vec, digits=5)), "\t") * "\t" *
                         join( string.(round.(clado_par_acceptfreq_vec ./clado_par_propfreq_vec, digits=5)), "\t")

        file_line=log_line_string

        log_file =  open(log_filename,"a")

        println(log_file, file_line)
        print(file_line)

        close(log_file)

   

        #setup proposal prob for proposing on a rate par of clado par

        clado_par_proposal_prob=1-rate_par_proposal_prob
    

        for iter in 1:iters

            #print("\n")
            #print(iter)
            #print("\n")
           


            par2upd=sample([1,2], Weights([rate_par_proposal_prob, clado_par_proposal_prob,]))
            #make the first iteration update rate par so Q_c gets made in loop


            if par2upd==1

                rate_pars_c ,
                edge_probs  , 
                post_clado_probs,
                rate_priorlL_c   ,
                prunelL_c   , 
                lL_c        ,
                rate_par_acceptfreq_vec ,
                rate_par_propfreq_vec =rate_par_upd(rate_pars_c      ,
                                                    tip_probs        ,
                                                    clado_Pmat_c       ,
                                                    lL_c             ,
                                                    prunelL_c        ,
                                                    clado_priorlL_c  ,
                                                    rate_priorlL_c   ,  
                                                    edge_probs         ,
                                                    post_clado_probs   ,                                                
                                                    Q_par_matrix     ,
                                                    Q_zeros          ,
                                                    Q_index_vec  ,
                                                    rate_prior_dists       ,
                                                    brs              ,
                                                    node_path        ,
                                                    rate_tuning_par_vec   , rate_par_acceptfreq_vec, rate_par_propfreq_vec ,
                                                    prior_only
                                                    )

            else 


            clado_probs_c ,
            edge_probs         ,
            post_clado_probs   ,
            clado_priorlL_c    ,
            prunelL_c          , 
            lL_c               ,
            clado_par_acceptfreq_vec ,
            clado_par_propfreq_vec =clado_par_upd(
                                              rate_pars_c      ,
                                              clado_probs_c     , 
                                              clado_Pmat_c           ,
                                              lL_c                   ,
                                              prunelL_c              ,
                                              clado_priorlL_c        ,
                                              rate_priorlL_c         ,  
                                              edge_probs         ,
                                              post_clado_probs   ,                                               
                                              clado_prior_dists       , 
                                              split_index_vec              ,
                                              sub_index_vec                ,
                                              Q_par_matrix     ,
                                              Q_zeros          ,
                                              Q_index_vec  ,
                                              tip_probs              ,
                                              brs                    ,
                                              node_path              ,
                                              clado_tuning_par_vec   ,
                                              clado_par_acceptfreq_vec     ,
                                              clado_par_propfreq_vec       ,
                                              prior_only
                                              )
            

            end
                                            
            if mod(iter, anc_state_sampling)==0

                
                anc_states, anc_clado_events=sample_anc_states(edge_probs,  post_clado_probs, root_probs, tip_probs, Q_zeros, rate_pars_c, Q_index_vec, Q_par_matrix, clado_Pmat_c, brs, node_path)

                #print(root_probs)
                #print("\n") 

                #print("\n")

                push!(anc_states_cache, join(string.(anc_states),"\t"))
                push!(anc_clado_cache, join(string.(anc_clado_events),"\t"))


                if mod(length(anc_states_cache), write_interval)==0

                    anc_states_log_file =  open(anc_states_log_filename,"a")

                    [println( anc_states_log_file, anc_states_cache[i]) for i in eachindex(anc_states_cache)]
                   
                    close( anc_states_log_file )

                    anc_states_cache=[]
                   
                end


                if mod(length(anc_clado_cache), write_interval)==0

                    anc_clado_log_file =  open( anc_clado_log_filename,"a")

                    [println( anc_clado_log_file, anc_clado_cache[i]) for i in eachindex(anc_clado_cache)]
                   
                    close( anc_clado_log_file )

                    anc_clado_cache=[]

                end

                
                #[println(log_file, chain_cache[i]) for i in eachindex(chain_cache)]



                #println(anc_clado_log_file, join(string.(anc_clado_events),"\t"))

                #println(anc_state_log_file, join(string.(anc_states),"\t"))
               # print("anc_states logged")
            end


            if mod(iter, iter_trims)==0
                file_line= string(iter) * "\t" *
                           string(round(lL_c, digits=4)) * "\t" * 
                           string(round(prunelL_c,digits=4)) * "\t" * 
                           string(round(rate_priorlL_c,digits=4)) * "\t" * 
                           string(round(clado_priorlL_c, digits=4)) * "\t" *
                           join( string.(round.(rate_pars_c, digits=5)), "\t") * "\t"  * 
                           join( string.(round.(clado_probs_c, digits=5)), "\t") * "\t"  * 
                           #string(round(clado_probs_c, digits=3)) * "\t"  * 
                           join( string.(round.(rate_par_acceptfreq_vec ./rate_par_propfreq_vec, digits=3)), "\t") * "\t" *
                           join( string.(round.(clado_par_acceptfreq_vec ./clado_par_propfreq_vec, digits=3)), "\t") 

                
                push!(chain_cache, file_line)


                if mod(length(chain_cache), write_interval)==0

                    log_file =  open(log_filename,"a")

                    [println(log_file, chain_cache[i]) for i in eachindex(chain_cache)]
                   
                    close(log_file)
                   
                    chain_cache=[]
                end
                
                #println(log_file, file_line)
                #print(file_line)
            end

        end 

    end   

    #end
#end
#end




function clado_multi_move(par::Float64, #parameter
    d::Float64  #tuning par
    )
    m = exp(d * (rand(Uniform(0,1)) - 0.5))
    par_prop=par * m
    hr=log(m)   

    return(par_prop, hr)

end



function multi_move(par::Float64, #parameter
                    d::Float64  #tuning par
    )
    m = exp(d * (rand(Uniform(0,1)) - 0.5))
    par_prop=par * m
    hr=log(m)

    return(par_prop, hr)

end
#prop <- as.vector(t(old)) * m


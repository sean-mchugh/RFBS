{
  
  createLayoutMatrix <- function(Nrow, Ncol, blockRows, blockCols, fillOrder = "byrow") {
    # Number of plots in a block
    nPlotsInBlock <- Nrow * Ncol
    
    # Total number of plots
    totalPlots <- nPlotsInBlock * blockRows * blockCols
    
    # Create an empty matrix
    layoutMatrix <- matrix(0, nrow = Nrow * blockRows, ncol = Ncol * blockCols)
    
    # Function to calculate plot number
    calculatePlotNumber <- function(br, bc, r, c) {
      if (fillOrder == "byrow") {
        return(((br - 1) * blockCols * nPlotsInBlock) + 
                 ((bc - 1) * Ncol) + 
                 ((r - 1) * Ncol * blockCols) + c)
      } else { # fillOrder == "bycolumn"
        return(((bc - 1) * blockRows * nPlotsInBlock) + 
                 ((br - 1) * Nrow) + 
                 ((c - 1) * Nrow * blockRows) + r)
      }
    }
    
    # Fill the matrix
    for (blockRow in 1:blockRows) {
      for (blockCol in 1:blockCols) {
        for (r in 1:Nrow) {
          for (c in 1:Ncol) {
            layoutMatrix[Nrow * (blockRow - 1) + r, Ncol * (blockCol - 1) + c] <- 
              calculatePlotNumber(blockRow, blockCol, r, c)
          }
        }
      }
    }
    
    return(layoutMatrix)
  }  
  
  t_col <- function(color, percent = 50, name = NULL) {
    #      color = color name
    #    percent = % transparency
    #       name = an optional name for the color
    
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    ## Save the color
    invisible(t.col)
  }
  
  
  get_aff_ratios=function(anc_states_post_burn, states_mat, strict_fund=F ){
    
    
    states_list=unlist(lapply(1:nrow(states_mat), function(i) paste(states_mat[i,],collapse = " ")))
    
    #made list objects to fill
    node_real_aff_ratios=list()
    node_fund_aff_ratios=list()
    
    node_non_aff_ratios=list()
    
    corner_real_aff_ratios=list()
    corner_fund_aff_ratios=list()
    
    corner_non_aff_ratios=list()
    
    #get node length (note this includes corners and nodes so must divide)
    node_length=ncol(anc_states_post_burn)/2
    
    for (nd in 1:node_length){
      
      cn= nd+node_length
      
      for (node in c(nd,cn)){
        
        st_freqs=table(anc_states_post_burn[,node])
        
        total_freq=sum(st_freqs)
        
        st_vecs=do.call(rbind, lapply(states_list, function(st) as.numeric(strsplit(st, " ")[[1]])))
        
        nbiomes=ncol(st_vecs)
        
        non_biome_freqs=vector(mode = "numeric",length = nbiomes)
        fund_biome_freqs=vector(mode = "numeric",length = nbiomes)
        real_biome_freqs=vector(mode = "numeric",length = nbiomes)
        
        
        
        for (i in 1:length(st_freqs)){
          
          st=as.numeric(names(st_freqs[i]))
          non_biomes=(1:nbiomes)[st_vecs[st,]==0]
          
          # strict fund is good for finding the exact number of specific support for fundamental affinities, not strict is good if trying to plot affinities, as if it is a real affinity by default it has to have a fundamental as well
          if(strict_fund==T){
            
            fund_biomes=(1:nbiomes)[st_vecs[st,]==1]
            
          }else{
            fund_biomes=(1:nbiomes)[st_vecs[st,]>0]
          }
          
          real_biomes=(1:nbiomes)[st_vecs[st,]>1]
          
          
          non_biome_freqs[non_biomes]=non_biome_freqs[non_biomes]+st_freqs[i]
          
          fund_biome_freqs[fund_biomes]=fund_biome_freqs[fund_biomes]+st_freqs[i]
          real_biome_freqs[real_biomes]=real_biome_freqs[real_biomes]+st_freqs[i]
          
        }
        
        if(node==cn){
          
          corner_non_aff_ratios[[nd]]=non_biome_freqs/total_freq
          corner_fund_aff_ratios[[nd]]=fund_biome_freqs/total_freq
          corner_real_aff_ratios[[nd]]=real_biome_freqs/total_freq
          
          #print(fund_biome_freqs)
          
        }else {
          
          node_non_aff_ratios[[nd]]=non_biome_freqs/total_freq
          node_fund_aff_ratios[[nd]]=fund_biome_freqs/total_freq
          node_real_aff_ratios[[nd]]=real_biome_freqs/total_freq
          
        }  
      }
    }
    
    node_non_aff_ratios_mat  =do.call(rbind,node_non_aff_ratios  ) 
    
    node_fund_aff_ratios_mat  =do.call(rbind,node_fund_aff_ratios  ) 
    node_real_aff_ratios_mat  =do.call(rbind,node_real_aff_ratios  )
    corner_non_aff_ratios_mat=do.call(rbind,corner_non_aff_ratios)
    
    corner_fund_aff_ratios_mat=do.call(rbind,corner_fund_aff_ratios)
    corner_real_aff_ratios_mat=do.call(rbind,corner_real_aff_ratios)
    
    
    
    
    
    
    
    return(list(node_non=node_non_aff_ratios_mat, node_fund=node_fund_aff_ratios_mat, node_real=node_real_aff_ratios_mat,
                corner_non=corner_non_aff_ratios_mat, corner_fund=corner_fund_aff_ratios_mat, corner_real=corner_real_aff_ratios_mat))
    
  }
  
  get_clado_moves=function(tree, anc_states_post_burn){
    
    clado_moves=vector()
    
    for (i in 1:length(unique(tree$edge[,1]))){
      print(i)
      
      nd=i+length(tree$tip.label)  
      row_index <- which(tree$edge[,1]==nd, arr.ind = TRUE)
      child <- tree$edge[row_index,2]
      
      clado_moves=append(clado_moves,apply(cbind(anc_states_post_burn[,nd], anc_states_post_burn[,child+node_length]),1,paste,collapse=" "))
      
      
    }
    
    return(clado_moves)
    
  }
  
  clado2cladotype=function(clado_move,states_mat){
    
    
    clado_states=states_mat[as.numeric(str_split(clado_move," ")[[1]]),]
    matching_R=all(clado_states[1,]==clado_states[2,])
    matching_L=all(clado_states[1,]==clado_states[3,])
    total_match=matching_R+matching_L
    if(total_match==2){
      
      clado_move_type="full_inheritance"
      
      
    }else if(total_match==1){
      
      allo_R=sum(clado_states[2,clado_states[1,]!=clado_states[2,]]==1)/length(clado_states[2,clado_states[1,]!=clado_states[2,]]==1)
      allo_L=sum(clado_states[3,clado_states[1,]!=clado_states[3,]]==1)/length(clado_states[2,clado_states[1,]!=clado_states[3,]]==1)
      
      eco_R=sum(clado_states[2,clado_states[1,]!=clado_states[2,]]==0) /length(clado_states[2,clado_states[1,]!=clado_states[2,]]==0)
      eco_L=sum(clado_states[3,clado_states[1,]!=clado_states[3,]]==0) /length(clado_states[2,clado_states[1,]!=clado_states[3,]]==0)
      
      
      
      if (eco_R ==1 || eco_L ==1){
        
        clado_move_type="partial_eco"
        
      }else if(allo_R ==1 || allo_L ==1){
        clado_move_type="partial_allo"
        
        
      }
    }else if(total_match==0){
      
      allo_R=sum(clado_states[2,clado_states[1,]!=clado_states[2,]]==1)/length(clado_states[2,clado_states[1,]!=clado_states[2,]]==1)
      allo_L=sum(clado_states[3,clado_states[1,]!=clado_states[3,]]==1)/length(clado_states[2,clado_states[1,]!=clado_states[3,]]==1)
      
      eco_R=sum(clado_states[2,clado_states[1,]!=clado_states[2,]]==0) /length(clado_states[2,clado_states[1,]!=clado_states[2,]]==0)
      eco_L=sum(clado_states[3,clado_states[1,]!=clado_states[3,]]==0) /length(clado_states[2,clado_states[1,]!=clado_states[3,]]==0)
      
      
      
      if (eco_R ==1 && eco_L ==1){
        
        clado_move_type="full_eco"
        
      }else if(allo_R ==1 && allo_L ==1){
        clado_move_type="full_allo"
        
        
      }else{
        
        return("error, didnt fill any category")
      }
      
    }
    
    print("type found")  
    
    return(clado_move_type)
  }
  
  recurse_gsub=function( patterns, replace="", strings){
    
    for( i in 1:length(patterns)){
      
      strings=gsub( patterns[[i]],replace, strings,fixed=FALSE,)
    }
    
    return(strings)
    
  }
  
  filelist2simstring=function(files){
    
    #files=subdir_file_list$tip_states_sims
    #femove all file paths and leave just the lowest directory file name
    files=lapply(files, function(f) str_split(f, "/")[[1]][[ length(str_split(f, "/")[[1]]) ]])
    files
    
    string_bits=c("anc_states.txt",
                  "_anc_acc",
                  "_anc_aff_post_sup",
                  "RFBS_",
                  "DEC_",
                  "_true_anc_aff",
                  "anc_clados.txt",
                  "anc_states.txt",
                  "_post_pred",
                  "___R_tree_[0-9]+.txt",
                  "tip_states.txt"
                  
    )
    
    #,
    #               "anc_"
    #              ,"_post_sup",
    #               "_acc",
    #               ".txt",
    #               "states",
    #               "_true",
    #               "_aff",
    #               "_clados",
    #               "_states",
    #               "_anc",
    #               "DEC_",
    #               "RFBS_",
    #               "_log")
    
    
    sim_pars_strings = recurse_gsub(string_bits, "", unlist(files))
    
    sim_pars_strings 
    #cleans all info from file names except the run number and simulating parammeters (runnum__rate_pars__cladopars)
    #sim_pars_strings =   gsub("_post_sup","", 
    #                          gsub("_acc","",
    #                            gsub(".txt","",
    #                               gsub("states","",
    #                                    gsub("_true","",
    #                                         gsub("_aff","", 
    #                                              gsub("_clados","",
    #                                                   gsub("_states","",
    #                                                        gsub("anc_","", 
    #                                                             gsub("DEC_","",
    #                                                                  gsub("RFBS_","", 
    #                                                                       gsub("_log","", unlist(files),fixed=FALSE),fixed=FALSE)))))))))))
    #
    #
    sim_pars_strings_= lapply(sim_pars_strings, function(f) unlist(str_split(f, "__"))  )
    
    #some labels have duplicates of run number, clean those out for the time being
    for (sim in 1:length(sim_pars_strings_)){
      
      sim_pars_strings_[[sim]][[1]]= unique(str_split_1(sim_pars_strings_[[sim]][[1]], "_"))
      
      sim_pars_strings[[sim]]=paste( unlist(sim_pars_strings_[[sim]]), collapse="__")
      
    }
    
    return(sim_pars_strings)
  }
  
  check_RFBSDEC_simfile_counts=function(sim_pars_strings){
    
    
    bad_runs=list()
    
    sim_pars_strings_unique=sort(unique(sim_pars_strings))
    
    dup_check_vec=unlist(lapply(sim_pars_strings_unique, function(sim) sum(sim==sim_pars_strings)))
    
    for(sim in 1:length(sim_pars_strings_unique)){
      
      if(dup_check_vec[[sim]]< max(dup_check_vec) ){
        
        print(paste("missing run for sim ",sim_pars_strings_unique[[sim]] ))
        
        bad_runs[[length(bad_runs)+1]]=as.numeric(str_split_1(sim_pars_strings_unique[[sim]], "_")[[1]])
        
        
      }
      
      
    }
    
    return(unlist(bad_runs))
    
  }
  
  clean_files=function(file_list){
    
    if (!is.list(file_list)){
      files_list=list(file_list)
    }
    
    bad_sims=unique(unlist(lapply(file_list, function(files) check_RFBSDEC_simfile_counts( filelist2simstring(c(files)) ) )))
    
    
    clean_files=  lapply(file_list, function(files)  files[!as.numeric(unlist(lapply(1:length(files), function(f) str_split_1( str_split_1( files[[f]], "/")[[length(str_split_1(files[[f]], "/"))]], "_" )[[1]] ))) %in% bad_sims] )
    
    return(clean_files)
    
  }
  
  sort_files=function(files, sep_string, files_name=""){
    
    
    n_sep=str_count(files[[1]],pattern = "/")
    file_mat=str_split_fixed(files, "/", n=n_sep)
    
    #matched_files=lapply(sep_string, function(string) grep(string,file_mat[,3]))
    
    #print(paste("found files for:", sep_string[unlist(lapply(matched_files, length))>0] ))
    
    #sort(unlist(  matched_files ))
    
    sep_file_mat_list = lapply(1:length(sep_string), function(i) file_mat[grep(sep_string[[i]],file_mat[,n_sep]),] )
    
    
    sep_file_mat_list = lapply(1:length(sep_string), function(s) as.vector(do.call(rbind, lapply(1:nrow(sep_file_mat_list[[s]]), function(i) paste(sep_file_mat_list[[s]][i,], collapse="/")))))
    
    names(sep_file_mat_list)=paste(files_name, sep_string, sep="")
    
    return(sep_file_mat_list )
    
    
  }
  
  getsubdir_files=function(dir){
    
    #pull sub directories paths
    subdir_list  = list.dirs(dir, recursive =F)
    #pull subdirfectory names
    subdir_names = list.dirs(dir, recursive =F,full.names = F)
    #make empty list to fille with subdirectory files
    subdir_file_list= list()
    
    for (d in 1:length(subdir_list)){
      subdir_file_list[[d]] = list.files(subdir_list[[d]], full.names=T)
    }
    
    #file groups correspond to subdirectory names
    names(subdir_file_list)=subdir_names
    
    #clean in case any list is missing a simulation run, and then drop all file with that run
    #test
    #subdir_file_list$anc_acc=subdir_file_list$anc_acc[-80]
    
    cleansubdir_files_list = clean_files(subdir_file_list)
    
    
    sep_string=c("RFBS", "DEC")
    
    
    
    
    sort_files(cleansubdir_files_list[[1]],   sep_string)
    
    
    
    return(cleansubdir_files_list)
    
  }
  
  aff_list2df=function(aff_ratios){
    node_aff_ratios_df=do.call(cbind, lapply(1:3, function(b)   cbind(aff_ratios$node_non[,b],   aff_ratios$node_fund[,b],    aff_ratios$node_real[,b])))
    corner_aff_ratios_df=do.call(cbind, lapply(1:3, function(b) cbind(aff_ratios$corner_non[,b], aff_ratios$corner_fund[,b],  aff_ratios$corner_real[,b])))
    
    full_aff_ratios=rbind(node_aff_ratios_df, corner_aff_ratios_df)
    
    
    
    return(full_aff_ratios)
    
    
    
    
  }
  
  get_affdf_colind=function(aff_df){
    
    nbiomes= ncol(aff_df)/3
    ntips=(nrow(aff_df)+2)/4
    deadrow=3*ntips
    nodes_ind= (ntips+1):((ntips*2)-1)

    col_indices=list(
      "non"= seq( from=1, by=3, length.out=nbiomes), # columns for non affinities, biome affinities are column bound, with the first col of each biome being the non affinity
      "fund"= seq( from=2, by=3, length.out=nbiomes), # second column of every biome is a fundamental affinity probability
      "real"= seq( from=3, by=3, length.out=nbiomes), # third column of every biome is the realized affinity probability 
      "unocc" =sort( c( seq( from=1, by=3, length.out=nbiomes),seq( from=2, by=3, length.out=nbiomes))),
      "aff"  =sort(  c(seq( from=2, by=3, length.out=nbiomes), seq( from=3, by=3, length.out=nbiomes)))
    )
    
    #calculate the dead row, corner affinities dont exist at the root 
    deadrow=ntips*3
    
    row_indices=list(
      "tips"=1:ntips,   #rows for tip affinities
      "nodes"=(ntips+1):((ntips*2)-1), #rows for node affinities 
      "corner"=((max(nodes_ind)+2):nrow(aff_df))[-deadrow] # rows for corner affintiies (post cladogenesis immediately after a node)
    )
    
    return( list("c"=col_indices, "r"= row_indices ))
    
  }
  
  plot_affinity_ratios=function(RFBS_aff_ratios, DEC_aff_ratios){
    
    DEC_all_node_0   = do.call(rbind, lapply(DEC_aff_ratios, function(sim) sim$node_non[-(1:ntips),]))
    DEC_all_node_2   = do.call(rbind, lapply(DEC_aff_ratios, function(sim) sim$node_real[-(1:ntips),]))
    DEC_all_corner_0 = do.call(rbind, lapply(DEC_aff_ratios, function(sim) sim$corner_non[-(1:ntips),]))
    DEC_all_corner_2 = do.call(rbind, lapply(DEC_aff_ratios, function(sim) sim$corner_real[-(1:ntips),]))
    
    RFBS_all_node_0   = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$node_non[-(1:ntips),]))
    RFBS_all_node_1   = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$node_fund[-(1:ntips),]))
    RFBS_all_node_2   = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$node_real[-(1:ntips),]))
    RFBS_all_corner_0 = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$corner_non[-(1:ntips),]))
    RFBS_all_corner_1 = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$corner_fund[-(1:ntips),]))
    RFBS_all_corner_2 = do.call(rbind, lapply(RFBS_aff_ratios, function(sim) sim$corner_real[-(1:ntips),]))
    
    
    {
      biome_labels=c("tropical", "warm temp", "cold temp")
      
      par(mfrow=c(3,2))
      
      for (b in 1:3){
        #boxplot(aff_ratios_list[[1]]$node_non[-(1:ntips),b],   aff_ratios_list[[2]]$node_non[-(1:ntips),b],names=c("DEC","RFBS" ), main=paste("non aff support",biome_labels[[b]]))
        #boxplot(aff_ratios_list[[1]]$node_real[-(1:ntips),b],   aff_ratios_list[[2]]$node_real[-(1:ntips),b], names=c("DEC","RFBS" ), main=paste("real aff support",biome_labels[[b]]))
        
        #vioplot(aff_ratios_list[[1]]$node_non[-(1:ntips),b],   aff_ratios_list[[2]]$node_non[-(1:ntips),b],names=c("DEC","RFBS" ), main=paste("non aff support",biome_labels[[b]]))
        #vioplot(aff_ratios_list[[1]]$node_real[-(1:ntips),b],   aff_ratios_list[[2]]$node_real[-(1:ntips),b], names=c("DEC","RFBS" ), main=paste("real aff support",biome_labels[[b]]))
        
        
        #compare exact affinities 
        vioplot( DEC_all_node_0[,b],                   0,  DEC_all_node_2[,b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "nodes"), side="left", col=2, colMed = 2, lineCol = 2)
        vioplot(RFBS_all_node_0[,b], RFBS_all_node_1[,b], RFBS_all_node_2[,b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=T, col=4, colMed = 4, lineCol = 4)
        
        vioplot( DEC_all_corner_0[,b],                     0,  DEC_all_corner_2[,b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "corners"), side="left",  col=2, colMed = 2, lineCol = 2)
        vioplot(RFBS_all_corner_0[,b], RFBS_all_corner_1[,b], RFBS_all_corner_2[,b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=T, col=4, colMed = 4, lineCol = 4)
        
        
        #compare support for any affinity at all across nodes
        #vioplot( DEC_aff_ratios[[i]]$node_non[-(1:ntips),b],    DEC_aff_ratios[[i]]$node_real[-(1:ntips),b], names=c("Non-affinity", "Affinity (Fund or Real)" ), main=paste(biome_labels[[b]], "nodes"), side="left", col=2, colMed = 2, lineCol = 2)
        #vioplot(RFBS_aff_ratios[[i]]$node_non[-(1:ntips),b],   RFBS_aff_ratios[[i]]$node_fund[-(1:ntips),b], names=c("Non-affinity", "Affinity (Fund or Real)" ), side="right",add=T, col=4, colMed = 4, lineCol = 4)
        #
        #vioplot(DEC_aff_ratios[[i]]$corner_non[-(1:ntips),b],     DEC_aff_ratios[[i]]$corner_real[-(1:ntips),b], names=c("Non-affinity","Affinity (Fund or Real)" ), main=paste(biome_labels[[b]], "corners"), side="left",  col=2, colMed = 2, lineCol = 2)
        #vioplot(RFBS_aff_ratios[[i]]$corner_non[-(1:ntips),b],   RFBS_aff_ratios[[i]]$corner_fund[-(1:ntips),b], names=c("Non-affinity","Affinity (Fund or Real)" ), side="right",add=T, col=4, colMed = 4, lineCol = 4)
        
        
        
        #vioplot( DEC_all_node_0[,b],    DEC_all_node_2[,b], names=aff_names, main=paste(biome_labels[[b]], "nodes"), side="left", col=2, colMed = 2, lineCol = 2)
        #vioplot(RFBS_all_node_0[,b],    RFBS_all_node_2[,b], names=aff_names, side="right",add=T, col=4, colMed = 4, lineCol = 4)
        #
        #vioplot( DEC_all_corner_0[,b],    DEC_all_corner_2[,b], names=aff_names, main=paste(biome_labels[[b]], "corners"), side="left", col=2, colMed = 2, lineCol = 2)
        #vioplot(RFBS_all_corner_0[,b],    RFBS_all_corner_2[,b], names=aff_names, side="right",add=T, col=4, colMed = 4, lineCol = 4)
        
        
        #boxplot(aff_ratios_list[[1]]$node_non[-(1:ntips),b],   0, aff_ratios_list[[1]]$node_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "nodes"), side="left", col=2, colMed = 2, lineCol = 2)
        #boxplot(aff_ratios_list[[2]]$node_non[-(1:ntips),b],   aff_ratios_list[[2]]$node_fund[-(1:ntips),b], aff_ratios_list[[2]]$node_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=F, col=4, colMed = 4, lineCol = 4)
        #boxplot(aff_ratios_list[[1]]$corner_non[-(1:ntips),b],   0, aff_ratios_list[[1]]$corner_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "corners"), side="left",  col=2, colMed = 2, lineCol = 2)
        #boxplot(aff_ratios_list[[2]]$corner_non[-(1:ntips),b],   aff_ratios_list[[2]]$corner_fund[-(1:ntips),b], aff_ratios_list[[2]]$corner_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=F, col=4, colMed = 4, lineCol = 4)
        
        
        #boxplot(aff_ratios_list[[1]]$node_non[-(1:ntips),b],   0, aff_ratios_list[[1]]$node_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "nodes"), side="left", col=2, colMed = 2, lineCol = 2)
        #boxplot(aff_ratios_list[[2]]$node_non[-(1:ntips),b],   aff_ratios_list[[2]]$node_fund[-(1:ntips),b], aff_ratios_list[[2]]$node_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=F, col=4, colMed = 4, lineCol = 4)
        #boxplot(aff_ratios_list[[1]]$corner_non[-(1:ntips),b],   0, aff_ratios_list[[1]]$corner_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), main=paste(biome_labels[[b]], "corners"), side="left",  col=2, colMed = 2, lineCol = 2)
        #boxplot(aff_ratios_list[[2]]$corner_non[-(1:ntips),b],   aff_ratios_list[[2]]$corner_fund[-(1:ntips),b], aff_ratios_list[[2]]$corner_real[-(1:ntips),b], names=c("Non-affinity","Fund-Affinity", "Real-Affinity" ), side="right",add=F, col=4, colMed = 4, lineCol = 4)
        
        
      }
      
      par(xpd=T)
      legend("bottom", inset=c(5, -.5), legend=c("DEC", "RFBS"), fill=c(2,4), cex = 1, horiz = T)
      
    }
    
    
  }
  
  tcol <- function(color, percent = 50, name = NULL) {
    #      color = color name
    #    percent = % transparency
    #       name = an optional name for the color
    
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    ## Save the color
    invisible(t.col)
  }
  
  
  get_aa_affs_from_dir=function(dir){
    
    dirfile_list   = getsubdir_files(dir)
    DEC_RFBS_logfiles = sort_files(dirfile_list$logs, c("RFBS", "DEC"), "logs")
    DEC_log_files  = sort_files(DEC_RFBS_logfiles$logsDEC, c("log$", "log_anc_states", "log_anc_clados"))
    RFBS_log_files = sort_files(DEC_RFBS_logfiles$logsRFBS, c("log$", "log_anc_states", "log_anc_clados"))
    
    #RFBS_as_post_list=lapply(RFBS_log_files$log_anc_states, function(sim) read.table(sim, header=F))
    
    dirfile_list$anc_acc 
    
    DEC_RFBS_aafiles = sort_files(dirfile_list$anc_acc , c("RFBS", "DEC"), "aa")
    DEC_aafiles  = sort_files(DEC_RFBS_aafiles$aaDEC, c("anc_aff_post_sup", "anc_acc"))
    RFBS_aafiles = sort_files(DEC_RFBS_aafiles$aaRFBS, c("anc_aff_post_sup", "anc_acc"))
    
    
    RFBS_states_space= read.table("rf_states.txt")
    DEC_states_space=read.table("DEC_states.txt")
    
    
    
    RFBS_aff_ratios_df=list()
    DEC_aff_ratios_df=list()
    true_aff_df=list()
    true_states=list()
    true_aff_list=list()
    
    for(i in 1:length(DEC_log_files$log_anc_states)){
      
      print(i)
      RFBS_aff_ratios_df[[i]]=read.table(RFBS_aafiles$anc_aff_post_sup[[i]])
      
      DEC_aff_ratios_df[[i]]=read.table(DEC_aafiles$anc_aff_post_sup[[i]])
      
      
      true_states[[i]]  =  na.omit(as.numeric(str_split_1(read.table(dirfile_list$anc_states_sims[[i]])[[1]], "\t")))
      true_aff_list[[i]]  =  get_aff_ratios(as.data.frame(t( true_states[[i]]), byrow = T) , RFBS_states_space, strict_fund=T)
      true_aff_df[[i]]  =  aff_list2df(true_aff_list[[i]])
    }
    
    
    
    return(list("RFBS_aff_post"=RFBS_aff_ratios_df,
                "DEC_aff_post"=  DEC_aff_ratios_df,
                "true_aff" =true_aff_df,
                "true_states"  = true_states
    ))
    
  }
  
}


library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(phytools)
library(phylobase) #to get node coords
library(plotrix)
library(vioplot)

{
  setwd("/Volumes/michael.landis/Active/RFBS_RIS")
  
  workdir=getwd()
  
  corner_coord=list.files(system.file("extdata/a_scripts", package="BioGeoBEARS"), full.names=TRUE)
  
  
  #dir_name="bvib_0301_runs_GOOD/Bvib_3nB_LN_0.1_1.0_1000000_eco_allo_clado_2g_1g_1l_rf_gl_ds"
  dir_names=c("BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.25]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.5]_rf_gl_ds",
              #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.75]_rf_gl_ds",
              "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds")
  
  
  dir_names=c("saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_1.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_1.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_1.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_1.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_5.0_r12_5.0_r02_0.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_5.0_r12_5.0_r02_0.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_5.0_r12_5.0_r02_5.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_1.0_r21_5.0_r12_5.0_r02_5.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_5.0_r21_1.0_r12_5.0_r02_0.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_5.0_r21_1.0_r12_5.0_r02_0.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_5.0_r21_1.0_r12_5.0_r02_5.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_1.0_r01_5.0_r21_1.0_r12_5.0_r02_5.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_1.0_r21_5.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_1.0_r21_5.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_1.0_r21_5.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_1.0_r21_5.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_5.0_r21_1.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_5.0_r21_1.0_r12_1.0_r02_0.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_5.0_r21_1.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_rf_gl_ds",
              "saved_BSim_runs/factorial/BSim_RFBS_DEC_comp_fixed150t_3nB_Exp1p0iter_300000_clado_Rtre_r10_5.0_r01_5.0_r21_1.0_r12_1.0_r02_1.0_splitP0.5_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds")
  
  
  label=c("all low (1.0)"                 ,"all low (1.0) uncertain tips", 
          "real high (5.0) fund low (1.0)","real high (5.0) fund low (1.0) uncertain tips",
          "gain high (5.0) loss low (1.0)","gain high (5.0) loss low (1.0) uncertain tips",
          "fund high (5.0) real low (1.0)","fund high (5.0) real low (1.0) uncertain tips",
          "loss high (5.0) gain low (1.0)","loss high (5.0) gain low (1.0) uncertain tips"
          
  )
  
  dir_ind=c(3, 4,
            7, 8,
            11, 12,
            15,16,
            19, 20)
  
  
  

  #list.dirs("saved_BSim_runs/BSim_runs/", recursive = F)
  
  
 dir_names=c("BSim_RFBS_DEC_comp50t_3nB_Exp0p5iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp50t_3nB_Exp1p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp50t_3nB_Exp5p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp150t_3nB_Exp5p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp500t_3nB_Exp0p5iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp500t_3nB_Exp1p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp500t_3nB_Exp5p0iter_500000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
 

 label=c("50 Tips \n 0.5 rate prior /",
         "50 Tips  \n 1.0 rate prior /",
         "50 Tips  \n 5.0 rate prior /",
         "150 Tips \n 0.5 rate prior /",
         "150 Tips \n 1.0 rate prior /",
         "150 Tips \n 5.0 rate prior /",
         "500 Tips \n 0.5 rate prior /",
         "500 Tips \n 1.0 rate prior /",
         "500 Tips \n 5.0 rate prior /")
 
 

 
 dir_names=c("BSim_RFBS_DEC_comp50t_3nB_Exp0p5iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             #"BSim_RFBS_DEC_comp50t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
             "BSim_RFBS_DEC_comp50t_3nB_Exp2p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
           # "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
           # "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
           # "BSim_RFBS_DEC_comp150t_3nB_Exp2p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
            "BSim_RFBS_DEC_comp500t_3nB_Exp0p5iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
            #"BSim_RFBS_DEC_comp500t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
            "BSim_RFBS_DEC_comp500t_3nB_Exp2p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
 
 
 label=c( "50 Tips \n 0.5 rate prior",
         #"50 Tips  \n 1.0 rate prior",
         "50 Tips  \n 2.0 rate prior",
         #"150 Tips \n 0.5 rate prior",
         #"150 Tips \n 1.0 rate prior",
         #"150 Tips \n 2.0 rate prior",
         "500 Tips \n 0.5 rate prior",
         #"500 Tips \n 1.0 rate prior",
         "500 Tips \n 2.0 rate prior")
 
 
 #dir=         "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds"
 #
 #dir_names=c("BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
 #
 #dir=dir_names[[1]]
 #
 #
 #dir_names=c(   "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.25]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.5]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][0.75]_rf_gl_ds",
 #               "saved_BSim_runs/uncertain_runs/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds")
 #
 #dir_names=c(   "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
 #               "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds")
 ##
 ##
 ##dir_ind=(c(1:9, 10))
 ##
 ##label=c("No ambiguous states \n in 100% tips", 
 ##        "1/3 ambiguous states \n in 25% tips",
 ##        "1/3 ambiguous states \n in 50% tips",
 ##        "1/3 ambiguous states \n in 75% tips",
 ##        "1/3 ambiguous states \n in 100% tips",
 ##        "2/3 ambiguous states \n in 25% tips",
 ##        "2/3 ambiguous states \n in 50% tips",
 ##        "2/3 ambiguous states \n in 75% tips",
 ##        "2/3 ambiguous states \n in 100% tips",
 ##        "all ambiguous states \n in 100% tips")
 ##
 #label=c("No ambiguous biomes \n in 100% tips","all ambiguous biomes \n in 100% tips",
 #        "-2 ambiguous biomes \n in 25% tips",  "-1 ambiguous state  \n in 25% tips",
 #        "-2 ambiguous biomes \n in 50% tips",  "-1 ambiguous state  \n in 50% tips",
 #        "-2 ambiguous biomes \n in 75% tips",  "-1 ambiguous state  \n in 75% tips",
 #        "-2 ambiguous biomes \n in 100% tips", "-1 ambiguous state  \n in 100% tips")
        
         
       
# dir_names=c(   "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
#                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
#               # "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.5]_rf_gl_ds",
#                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
#               # "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][1.0]_rf_gl_ds",
#                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
#               # "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.5]_rf_gl_ds",
#                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds",
#               # "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][1.0]_rf_gl_ds",
#                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][1.0][1.0]_rf_gl_ds")
# 
 
 # dir_names=c("saved_BSim_runs/Draft_2/BSim_RFBS_DEC_compPO_50t_3nB_Exp0p5iter_500000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")           
 
# dir_labels=c("No Missing Data",
#              "-2 ambiguous states of a maximum of 3 in 25% tips",
#              #"-2 ambiguous states of a maximum of 3 in 50% tips",
#              "-2 Ambiguous States in 75% Tips",
#             #   "-2 ambiguous states of a maximum of 3 in 100% tips",
#              "-1 ambiguous state of a maximum of 3 in 25% tips",
#             # "-1 ambiguous state of a maximum of 3 in 50% tips",
#              "-1 Ambiguous State in 75% Tips",
#              #"-1 ambiguous state of a maximum of 3 in 100% tips",
#              "All Ambiguous States in 100% Tips")#
 
 
         
         
  
  
}

burnin=0.5
RFBS_states_space= read.table("rf_states.txt")
DEC_states_space=read.table("DEC_states.txt")


get_aa_affs_from_dir=function(dir){

dirfile_list   = getsubdir_files(dir)
DEC_RFBS_logfiles = sort_files(dirfile_list$logs, c("RFBS", "DEC"), "logs")
DEC_log_files  = sort_files(DEC_RFBS_logfiles$logsDEC, c("log$", "log_anc_states", "log_anc_clados"))
RFBS_log_files = sort_files(DEC_RFBS_logfiles$logsRFBS, c("log$", "log_anc_states", "log_anc_clados"))

#RFBS_as_post_list=lapply(RFBS_log_files$log_anc_states, function(sim) read.table(sim, header=F))

dirfile_list$anc_acc 

DEC_RFBS_aafiles = sort_files(dirfile_list$anc_acc , c("RFBS", "DEC"), "aa")
DEC_aafiles  = sort_files(DEC_RFBS_aafiles$aaDEC, c("anc_aff_post_sup", "anc_acc"))
RFBS_aafiles = sort_files(DEC_RFBS_aafiles$aaRFBS, c("anc_aff_post_sup", "anc_acc"))





RFBS_aff_ratios_df=list()
DEC_aff_ratios_df=list()
true_aff_df=list()
true_states=list()
true_aff_list=list()

for(i in 1:length(DEC_log_files$log_anc_states)){
  
  print(i)
  RFBS_aff_ratios_df[[i]]=read.table(RFBS_aafiles$anc_aff_post_sup[[i]])
  
  DEC_aff_ratios_df[[i]]=read.table(DEC_aafiles$anc_aff_post_sup[[i]])
  
  
  true_states[[i]]  =  na.omit(as.numeric(str_split_1(read.table(dirfile_list$anc_states_sims[[i]])[[1]], "\t")))
  true_aff_list[[i]]  =  get_aff_ratios(as.data.frame(t( true_states[[i]]), byrow = T) , RFBS_states_space, strict_fund=T)
  true_aff_df[[i]]  =  aff_list2df(true_aff_list[[i]])
}



return(list("RFBS_aff_post"=RFBS_aff_ratios_df,
            "DEC_aff_post"=  DEC_aff_ratios_df,
            "true_aff" =true_aff_df,
            "true_states"  = true_states
            ))

}


affs_list=list()

for(d in 1:length(dir_names)){
  affs_list[[d]]=get_aa_affs_from_dir(dir_names[[d]])
}

#saveRDS(affs_list, "factorial_affinities.rds")
#affs_list=readRDS( "factorial_affinities.RDS")

getaff_stats=function(RFBS_aff_ratios_df,DEC_aff_ratios_df, true_aff_df, RFBS_states_space, DEC_states_space   ){
  #d=3
  #dir_names[[d]]
  
  #RFBS_aff_ratios_df = affs_list[[d]]$RFBS_aff_post
  #DEC_aff_ratios_df  = affs_list[[d]]$DEC_aff_post
  #true_aff_df        = affs_list[[d]]$true_aff
  
  RFBS2rprob=sum(RFBS_states_space[,1]==2) / nrow(RFBS_states_space)
  RFBS1rprob=sum(RFBS_states_space[,1]==1) / nrow(RFBS_states_space)
  RFBS0rprob=sum(RFBS_states_space[,1]==0) / nrow(RFBS_states_space)
  
  DEC2rprob=sum( DEC_states_space[,1]==2) / nrow(DEC_states_space)
  DEC0rprob =sum( DEC_states_space[,1]==0) / nrow(DEC_states_space)
  
  
  
  {
  RFBS_sim_acc=list()
  DEC_sim_acc=list()
  
  RFBS02_sim_acc=list()
  DEC02_sim_acc=list()
  
  
  
  RFBS2sim_acc=list()
  DEC2sim_acc=list()
  RFBS0sim_acc=list()
  DEC0sim_acc=list()
  RFBS1sim_acc=list()
  RFBS_DEC2=list()
  RFBS_DEC0=list()
  
  
  median_RFBS2sim_acc=list()
  median_DEC2sim_acc=list()
  median_RFBS0sim_acc=list()
  median_DEC0sim_acc=list()
  median_RFBS1sim_acc=list()
  
  
  rRFBS2_acc=list()
  rDEC2_acc=list()
  rRFBS0_acc=list()
  rDEC0_acc=list()
  rRFBS1_acc=list()
  rRFBS_DEC2=list()
  rRFBS_DEC0=list()
  RFBS2_acc=list()
  DEC2_acc=list()
  RFBS0_acc=list()
  DEC0_acc=list()
  RFBS1_acc=list()
  med_RFBS2_acc=list()
  med_DEC2_acc=list()
  med_RFBS0_acc=list()
  med_DEC0_acc=list()
  med_RFBS1_acc=list()
  
  
  RFBS_DEC2=list()
  RFBS_DEC0=list()
  
  
  RF_ind=get_affdf_colind(DEC_aff_ratios_df[[1]])
  
  RFrows=c(RF_ind$r$corner)
  
  RFBS2sim_meanprob=list()
  DEC2sim_meanprob =list()
  RFBS0sim_meanprob=list()
  DEC0sim_meanprob =list()
  
  rRFBSnonocc_acc=list()
  rDECnonocc_acc=list() 
  rRFBSaff_acc=list()   
  rDECaff_acc=list()    
  
  med_DECnonocc_acc=list() 
  med_RFBSnonocc_acc=list()
  med_DECaff_acc=list()    
  med_RFBSaff_acc=list() 
  
  RFBSnonocc_acc=list()    
  RFBSaff_acc=list()       
  DECnonocc_acc=list()     
  DECaff_acc=list()   
  
  median_RFBSnonoccsim_acc=list()
   median_DECnonoccsim_acc=list()
     median_RFBSaffsim_acc=list()
      median_DECaffsim_acc=list()
      RFBSnonoccsim_acc=list()
      DECnonoccsim_acc=list()
      RFBSaffsim_acc=list()
      DECaffsim_acc=list()
      
  
  for(i in 1:length(RFBS_aff_ratios_df)){
    
    print(i)
    #RFBS2_acc = RFBS_aff_ratios_df[[i]][RF_ind$r$nodes,RF_ind$c$real][true_aff_df[[i]][RF_ind$r$nodes,RF_ind$c$real]==1]/RFBS2rprob
    #DEC2_acc =  DEC_aff_ratios_df[[i]][RF_ind$r$nodes,RF_ind$c$real][true_aff_df[[i]][RF_ind$r$nodes,RF_ind$c$real]==1] / DEC2rprob
    #
    #RFBS2sim_acc[[i]]=(sum(RFBS2_acc)/length(RFBS2_acc))
    #DEC2sim_acc[[i]]= (sum(DEC2_acc)/length(DEC2_acc))
    #
    #RFBS0_acc=RFBS_aff_ratios_df[[i]][RF_ind$r$nodes,RF_ind$c$non][true_aff_df[[i]][RF_ind$r$nodes,RF_ind$c$non]==1]/RFBS0rprob
    #DEC0_acc = DEC_aff_ratios_df[[i]][RF_ind$r$nodes,RF_ind$c$non][true_aff_df[[i]][RF_ind$r$nodes,RF_ind$c$non]==1]/ DEC0rprob
    #
    #RFBS0sim_acc[[i]]= (sum(RFBS0_acc)/length(RFBS0_acc))
    #DEC0sim_acc[[i]] = (sum(DEC0_acc)/length(DEC0_acc))
    #
    #RFBS1_acc=RFBS_aff_ratios_df[[i]][RF_ind$r$nodes,RF_ind$c$fund][true_aff_df[[i]][RF_ind$r$nodes,RF_ind$c$fund]==1]
    #RFBS1sim_acc[[i]]= sum(RFBS1_acc)/length(RFBS1_acc)
    
    #RFBS02_sim_acc     = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,c(RF_ind$c$real,RF_ind$c$non) ]==1]#/RFBS2rprob
    #DEC02_sim_acc      = DEC_aff_ratios_df[[i]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,c(RF_ind$c$real,RF_ind$c$non) ]==1]#/RFBS2rprob
    RFBS2_acc[[i]]         = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$real]==1]                         #/RFBS2rprob
    DEC2_acc[[i]]          =  DEC_aff_ratios_df[[i]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$real]==1]                          #/ DEC2rprob
    
    RFBS0_acc[[i]]         = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$non][true_aff_df[[i]][RFrows,RF_ind$c$non]==1]                             #/RFBS0rprob
    DEC0_acc[[i]]          =  DEC_aff_ratios_df[[i]][RFrows,RF_ind$c$non][true_aff_df[[i]][RFrows,RF_ind$c$non]==1]                             #/ DEC0rprob
    
    RFBS1_acc[[i]]         = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$fund][true_aff_df[[i]][RFrows,RF_ind$c$fund]==1]                           #/RFBS1rprob
    #RFBS_occ_acc[[i]]         = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$real]==1]                         #/RFBS2rprob
    RFBSnonocc_acc[[i]]      = RFBS_aff_ratios_df[[i]][ RFrows,RF_ind$c$fund ][ true_aff_df[[i]][RFrows,RF_ind$c$non]==1 |true_aff_df[[i]][RFrows,RF_ind$c$fund]==1  ] + RFBS_aff_ratios_df[[i]][ RFrows,RF_ind$c$non ][ true_aff_df[[i]][RFrows,RF_ind$c$non]==1 |true_aff_df[[i]][RFrows,RF_ind$c$fund]==1  ]      
      #/RFBS2rprob
    RFBSaff_acc[[i]]         = RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$fund][  true_aff_df[[i]][RFrows,RF_ind$c$real]==1|true_aff_df[[i]][RFrows,RF_ind$c$fund]==1     ] + RFBS_aff_ratios_df[[i]][RFrows,RF_ind$c$real][  true_aff_df[[i]][RFrows,RF_ind$c$real]==1|true_aff_df[[i]][RFrows,RF_ind$c$fund]==1     ] 
    
    DECnonocc_acc[[i]]      =DEC_aff_ratios_df[[i]][ RFrows,RF_ind$c$non ][   true_aff_df[[i]][RFrows,RF_ind$c$non]==1 |true_aff_df[[i]][RFrows,RF_ind$c$fund]==1 ]                         #/RFBS2rprob
    DECaff_acc[[i]]         =DEC_aff_ratios_df[[i]][ RFrows,RF_ind$c$real ][  true_aff_df[[i]][RFrows,RF_ind$c$real]==1|true_aff_df[[i]][RFrows,RF_ind$c$fund]==1  ]                         #/RFBS2rprob
    

    
    med_RFBS0_acc[[i]]=median(RFBS0_acc[[i]])
    med_DEC0_acc[[i]] =median(DEC0_acc[[i]] )
    med_RFBS2_acc[[i]]=median(RFBS2_acc[[i]])
    med_DEC2_acc[[i]] =median(DEC2_acc[[i]] ) 
    med_RFBS1_acc[[i]]=median(RFBS1_acc[[i]])
    med_DECnonocc_acc[[i]] =median(DECnonocc_acc[[i]] )
    med_RFBSnonocc_acc[[i]]=median(RFBSnonocc_acc[[i]])
    med_DECaff_acc[[i]]    =median(DECaff_acc[[i]] ) 
    med_RFBSaff_acc[[i]]   =median(RFBSaff_acc[[i]])
    
    
    #rRFBS02_sim_acc[[i]]     =  lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,c(RF_ind$c$real,RF_ind$c$non) ]==1])#/RFBS2rprob
    #rDEC02_sim_acc[[i]]       =  lapply(1:length((RFBS_aff_ratios_df)), function(s)  DEC_aff_ratios_df[[s]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,c(RF_ind$c$real,RF_ind$c$non) ]==1])#/RFBS2rprob
    rRFBS2_acc[[i]]           =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$real]==1]                 )))        #/RFBS2rprob
    rDEC2_acc[[i]]            =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s)  DEC_aff_ratios_df[[s]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$real]==1]                 )))         #/ DEC2rprob
    rRFBS0_acc[[i]]           =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$non][true_aff_df[[i]][RFrows,RF_ind$c$non]==1]                   )))          #/RFBS0rprob
    rDEC0_acc[[i]]            =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s)  DEC_aff_ratios_df[[s]][RFrows,RF_ind$c$non][true_aff_df[[i]][RFrows,RF_ind$c$non]==1]                   )))          #/ DEC0rprob
    rRFBS1_acc[[i]]           =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$fund][true_aff_df[[i]][RFrows,RF_ind$c$fund]==1]                 )))          #/RFBS1rprob
    
    rRFBSnonocc_acc[[i]]        =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$unocc][true_aff_df[[i]][RFrows,RF_ind$c$unocc]==1]                   )))          #/RFBS0rprob
    rDECnonocc_acc[[i]]         =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s)  DEC_aff_ratios_df[[s]][RFrows,RF_ind$c$non][true_aff_df[[i]][RFrows,RF_ind$c$unocc]==1]                 )))         #/ DEC2rprob
    rRFBSaff_acc[[i]]           =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s) RFBS_aff_ratios_df[[s]][RFrows,RF_ind$c$aff][true_aff_df[[i]][RFrows,RF_ind$c$aff]==1]                 )))          #/RFBS1rprob
    rDECaff_acc[[i]]           =  mean(unlist(lapply(1:length((RFBS_aff_ratios_df)), function(s)  DEC_aff_ratios_df[[s]][RFrows,RF_ind$c$real][true_aff_df[[i]][RFrows,RF_ind$c$aff]==1]                   )))          #/ DEC0rprob
    
    
    #RFBS2sim_acc[[i]]=log(sum(RFBS2_acc)/length(RFBS2_acc))
    #DEC2sim_acc[[i]]= log(sum(DEC2_acc)/length(DEC2_acc))
    #RFBS_DEC2[[i]]=(log(sum(RFBS2_acc))-log(sum(DEC2_acc)))
    #
    #RFBS0sim_acc[[i]]= log(sum(RFBS0_acc)/length(RFBS0_acc))
    #DEC0sim_acc[[i]] = log(sum(DEC0_acc)/length(DEC0_acc))
    #RFBS_DEC0[[i]]=(log(sum(RFBS0_acc))-log(sum(DEC0_acc)))
    #
    #RFBS1sim_acc[[i]]= log(sum(RFBS1_acc))/length(RFBS1_acc)
    
    

    RFBS0sim_acc[[i]] =RFBS0_acc[[i]]/rRFBS0_acc[[i]]
    DEC0sim_acc[[i]] = DEC0_acc[[i]]/rDEC0_acc[[i]]
    RFBS2sim_acc[[i]]= RFBS2_acc[[i]]/rRFBS2_acc[[i]]
    DEC2sim_acc[[i]]=  DEC2_acc[[i]]/  rDEC2_acc[[i]]
    RFBS1sim_acc[[i]]= ((RFBS1_acc[[i]]/rRFBS1_acc[[i]]))
    
    median_RFBSnonoccsim_acc[[i]]= (RFBSnonocc_acc[[i]]/rRFBSnonocc_acc[[i]])
    median_DECnonoccsim_acc[[i]] = (DECnonocc_acc[[i]] /rDECnonocc_acc[[i]] )
    median_RFBSaffsim_acc[[i]]=    (RFBSaff_acc[[i]]   /rRFBSaff_acc[[i]]   )
    median_DECaffsim_acc[[i]]=     (DECaff_acc[[i]]    /rDECaff_acc[[i]]    )
    
    
    #RFBS_DEC2[[i]]=((sum(RFBS2_acc))-(sum(DEC2_acc)))
    
    #RFBS_DEC0[[i]]=((sum(RFBS0_acc[[i]]))-(sum(DEC0_acc[[i]])))
    median_RFBS2sim_acc[[i]]= median(RFBS2_acc[[i]]/rRFBS2_acc[[i]])
    median_DEC2sim_acc[[i]]=  median(DEC2_acc[[i]]/  rDEC2_acc[[i]])
    median_RFBS1sim_acc[[i]]= median(((RFBS1_acc[[i]]/rRFBS1_acc[[i]])))
    
    #RFBS_DEC2[[i]]=((sum(RFBS2_acc))-(sum(DEC2_acc)))
    
    median_RFBS0sim_acc[[i]]=  median(RFBS0_acc[[i]]/rRFBS0_acc[[i]])
    median_DEC0sim_acc[[i]] = median(DEC0_acc[[i]]/rDEC0_acc[[i]])
    
    median_RFBSnonoccsim_acc[[i]]= median(RFBSnonocc_acc[[i]]/rRFBSnonocc_acc[[i]])
    median_DECnonoccsim_acc[[i]] = median(DECnonocc_acc[[i]] /rDECnonocc_acc[[i]] )
    median_RFBSaffsim_acc[[i]]=    median(RFBSaff_acc[[i]]   /rRFBSaff_acc[[i]]   )
    median_DECaffsim_acc[[i]]=     median(DECaff_acc[[i]]    /rDECaff_acc[[i]]    )
    
    
    
    
  }
  
  return(list("aff_acc"=list("RFBS2_acc"=RFBS2_acc,
                             "DEC2_acc" =DEC2_acc ,  
                             "RFBS0_acc"=RFBS0_acc,
                             "DEC0_acc" =DEC0_acc ,
                             "RFBS1_acc"=RFBS1_acc,
                             "RFBSnonocc_acc"= RFBSnonocc_acc,
                             "DECnonocc_acc"= DECnonocc_acc,
                             "RFBSaff_acc"= RFBSaff_acc,
                             "DECaff_acc"= DECaff_acc),
        "random_aff_acc"=list("RFBS2_acc"=rRFBS2_acc,
                              "DEC2_acc" =rDEC2_acc ,
                              "RFBS0_acc"=rRFBS0_acc,
                              "DEC0_acc" =rDEC0_acc ,
                              "RFBS1_acc"=rRFBS1_acc,
                              "RFBSnonocc_acc" = rRFBSnonocc_acc,
                              "DECnonocc_acc"  = rDECnonocc_acc,
                              "RFBSaff_acc"    = rRFBSaff_acc,
                              "DECaff_acc"     = rDECaff_acc),
        
     "corrected_aff_acc"=list("RFBS2_acc"=RFBS2sim_acc,
                              "DEC2_acc" =DEC2sim_acc ,
                              "RFBS0_acc"=RFBS0sim_acc,
                              "DEC0_acc" =DEC0sim_acc ,
                              "RFBS1_acc"=RFBS1sim_acc,
                              "RFBSnonocc_acc" = RFBSnonoccsim_acc,
                              "DECnonocc_acc"  = DECnonoccsim_acc,
                              "RFBSaff_acc"    = RFBSaffsim_acc,
                              "DECaff_acc"     = DECaffsim_acc), 
  "median_aff_acc"=list("RFBS2_acc"= med_RFBS2_acc,
                    "DEC2_acc"     = med_DEC2_acc ,  
                    "RFBS0_acc"    = med_RFBS0_acc,
                    "DEC0_acc"     = med_DEC0_acc ,
                    "RFBS1_acc"    = med_RFBS1_acc,
                    "RFBSnonocc_acc" = med_RFBSnonocc_acc,
                    "DECnonocc_acc"  = med_DECnonocc_acc,
                    "RFBSaff_acc"    = med_RFBSaff_acc,
                    "DECaff_acc"     = med_DECaff_acc),
  "median_corrected_aff_acc"=list("RFBS2_acc"= median_RFBS2sim_acc,
                                  "DEC2_acc" = median_DEC2sim_acc ,
                                  "RFBS0_acc"= median_RFBS0sim_acc,
                                  "DEC0_acc" = median_DEC0sim_acc ,
                                  "RFBS1_acc"= median_RFBS1sim_acc,
                                  "RFBSnonocc_acc" = median_RFBSnonoccsim_acc,
                                  "DECnonocc_acc"  = median_DECnonoccsim_acc,
                                  "RFBSaff_acc"    = median_RFBSaffsim_acc,
                                  "DECaff_acc"     = median_DECaffsim_acc)
  ))

}

  
  
  #mean(unlist(RFBS2sim_acc))
#mean(unlist(RFBS1sim_acc))
#mean(unlist(RFBS0sim_acc))
#
#mean(unlist(rRFBS2_acc))
#mean(unlist(rRFBS1_acc))
#mean(unlist(rRFBS0_acc))
#
#
#mean(unlist(DEC2sim_acc))
#mean(unlist(DEC0sim_acc))
#
#density(unlist(RFBS02_sim_acc))
#
#RFBS2sim_accden=(density(unlist(RFBS2sim_acc)))
#RFBS1sim_accden=(density(unlist(RFBS1sim_acc)))
#RFBS0sim_accden=(density(unlist(RFBS0sim_acc)))
#DEC2sim_accden= (density(unlist(DEC2sim_acc )))
#DEC0sim_accden= (density(unlist(DEC0sim_acc )))
}


aff_acc_list=lapply(1:length(affs_list), function(d) getaff_stats(affs_list[[d]]$RFBS_aff_post,
             affs_list[[d]]$DEC_aff_post,
             affs_list[[d]]$true_aff, RFBS_states_space, DEC_states_space ))


d=3

#par(mfrow=c(1,3))

col_vec=c("dark blue", "light blue", "dark red", "pink")
col_vec=c("light blue", "pink")
border_vec= c("dark blue", "dark red")



dir_ind=1:10

dir_ind=(c(1,3,7,9))

dir_ind=1:4

label=c("50 Tips \n 0.5 rate prior",
        #"50 Tips  \n 1.0 rate prior /",
        "50 Tips  \n 2.0 rate prior",
        #"150 Tips \n 0.5 rate prior /",
        #"150 Tips \n 1.0 rate prior /",
        #"150 Tips \n 2.0 rate prior /",
        "500 Tips \n 0.5 rate prior",
        #"500 Tips \n 1.0 rate prior /",
        "500 Tips \n 2.0 rate prior")


#dir_ind=(c(1,4,8,10))
#
#label=c("no ambiguous state  \n in 100% tips", "-2 ambiguous states \n in 75% tips",
#        "-1 ambiguous state  \n in 75% tips","all ambiguous states \n in 100% tips")
#

#label=c("No Missing Data",
#             "-2 Ambiguous States in 25% tips",
#             #"-2 ambiguous states of a maximum of 3 in 50% tips",
#             "-2 Ambiguous States in 75% Tips",
#             #  "-2 ambiguous states of a maximum of 3 in 100% tips",
#             "-1 Ambiguous State  in 25% tips",
#             #"-1 ambiguous state of a maximum of 3 in 50% tips",
#             "-1 Ambiguous State in 75% Tips",
#             #"-1 ambiguous state of a maximum of 3 in 100% tips",
#             "All Ambiguous States in 100% Tips")
#
#label=c("No Ambiguous Data",
#        "Large Ambiguity Reduction in Some Tips",
#        #"-2 ambiguous states of a maximum of 3 in 50% tips",
#        "Large Ambiguity Reduction in Most Tips",
#        #  "-2 ambiguous states of a maximum of 3 in 100% tips",
#        "Small Ambiguity Reduction in Some Tips",
#        #"-1 ambiguous state of a maximum of 3 in 50% tips",
#        "Large Ambiguity Reduction in Most Tips",
#        #"-1 ambiguous state of a maximum of 3 in 100% tips",
#        "All Ambiguous Data in all Tips")
#
#
#
#dir_ind=c(1,3,2,5,4,6)
#######plot affs########

{

pdf(paste("3_RFnon_boxplots.pdf", sep="/"),width = 15,height = 15)  

par(mfrow=c(2,2), mai=rep(0.8,4))


for (i in 1:length(dir_ind)){
d=dir_ind[[i]]
  
boxplot(unlist(aff_acc_list[[d]]$median_aff_acc$RFBS2_acc), 
        #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc),
        unlist(aff_acc_list[[d]]$median_aff_acc$DEC2_acc ), 
      #  unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc ), 
        
        unlist(aff_acc_list[[d]]$median_aff_acc$RFBS1_acc),
        #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS1_acc),
        NA,
      #NA, 
        
        unlist(aff_acc_list[[d]]$median_aff_acc$RFBS0_acc) , 
        #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc), 
        unlist(aff_acc_list[[d]]$median_aff_acc$DEC0_acc )  ,
      #  unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc ), 
        
        #main=label[[i]] , #at=c(2,3, 4,5, 7,8,9,10,12,13,14, 15),#22,23), \
      
        at=c(2,3,
             5,6,
             8,9
        ),
        names = c("        Established \n       Affinity","",
                  "        Enabled \n       Affinity", "",
                  "        Non \n       Affinity", ""
                 
                  ), col=col_vec,border=border_vec, ylab="True Affinity Accuracy",
       pch=20, notch=F, xaxt="n",  ylim=c(0,1),
      cex.main=2 ,
      cex.lab=2)

       mtext(c("        Established \n       Affinity","",
               
               "        Enabled \n       Affinity", "",
               
               "        Non \n       Affinity", ""),
             1, 2, 
             at=c(2,3,
                  5,6,
                  8,9)
             )
       
       
       
       
       boxplot(unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc), 
               #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc),
                 unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc ), 
               #  unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc ), 
               
               unlist(aff_acc_list[[d]]$random_aff_acc$RFBS1_acc),
               #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS1_acc),
               NA,
               #NA, 
               
               unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc) , 
               #unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc), 
               unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc )  ,
               at=c(2,3,
                    5,6,
                    8,9
               ),
               add=T, col=t_col("gray", percent = 20),border = t_col("gray", percent = 20),outline=F, xaxt="n")
       
}

dev.off()


}

###plot aff/occ#################################################
{

pdf(paste("3_aff_occ_acc_boxplots.pdf", sep="/"),width = 15,height = 15)  

par(mfrow=c(2,2), mai=c(0.8,0.8,0.8,0.8))


for (i in 1:length(dir_ind)){
  d=dir_ind[[i]]
  
  boxplot(unlist(aff_acc_list[[d]]$median_aff_acc$RFBSaff_acc), 
         # unlist(aff_acc_list[[d]]$random_aff_acc$RFBSaff_acc),
            unlist(aff_acc_list[[d]]$median_aff_acc$DECaff_acc), 
          #  unlist(aff_acc_list[[d]]$random_aff_acc$DECaff_acc), 
          
          unlist(aff_acc_list[[d]]$median_aff_acc$RFBS0_acc),
         # unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc),
          unlist(aff_acc_list[[d]]$median_aff_acc$DEC0_acc),
          #unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc),
          
          
          unlist(aff_acc_list[[d]]$median_aff_acc$RFBS2_acc) , 
         # unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc), 
          unlist(aff_acc_list[[d]]$median_aff_acc$DEC2_acc)  ,
         # unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc), 
          
          unlist(aff_acc_list[[d]]$median_aff_acc$RFBSnonocc_acc) , 
         # unlist(aff_acc_list[[d]]$random_aff_acc$RFBSnonocc_acc), 
          unlist(aff_acc_list[[d]]$median_aff_acc$DECnonocc_acc)  ,
         # unlist(aff_acc_list[[d]]$random_aff_acc$DECnonocc_acc), 
          
          
          
          
          #main=label[[i]] , 
         at=c(2,3,
              5,6,
              8,9,
              11,12
               ),
          #22,23), 
          #names = c("Biome","Affinity",
          #          
          #          "Biome", "Non Affinity",
          #          
          #           "Occupied", "Biome",
          #          
          #           "unOccupied", "Biome"
          #          ), 
          col=col_vec,border=border_vec, ylab="True Affinity Accuracy",
          pch=20,
          ylim=c(0,1), xaxt="n",
         cex.main=2 ,
         cex.lab=2)
  
  mtext(c("        Biomel \n       Affinity","",
          
          "        Biome \n      Non-Affinity", "",
          
          "        Biome \n       Established", "",
          "        Biome \n       Unestablished", ""),
        1, 2, 
        at=c(2,3,
             5,6,
             8,9,
             11,12)
  )
  
  
  boxplot(unlist(aff_acc_list[[d]]$random_aff_acc$RFBSaff_acc), 
          # unlist(aff_acc_list[[d]]$random_aff_acc$RFBSaff_acc),
          unlist(aff_acc_list[[d]]$random_aff_acc$DECaff_acc), 
          #  unlist(aff_acc_list[[d]]$random_aff_acc$DECaff_acc), 
          
          unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc),
          # unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc),
          unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc),
          #unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc),
          
          
          unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc) , 
          # unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc), 
          unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc)  ,
          # unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc), 
          
          unlist(aff_acc_list[[d]]$random_aff_acc$RFBSnonocc_acc) , 
          # unlist(aff_acc_list[[d]]$random_aff_acc$RFBSnonocc_acc), 
          unlist(aff_acc_list[[d]]$random_aff_acc$DECnonocc_acc)  ,
          # unlist(aff_acc_list[[d]]$random_aff_acc$DECnonocc_acc), 
          
          
          
          
          #main=label[[i]] ,
          at=c(2,3,
                                 5,6,
                                 8,9,
                                 11,12
          ),
          #22,23), 
          #names = c("Biome","Affinity",
          #          
          #          "Biome", "Non Affinity",
          #          
          #           "Occupied", "Biome",
          #          
          #           "unOccupied", "Biome"
          #          ), 
          col=t_col("gray", percent = 20),border=t_col("gray", percent = 20), ylab="True Affinity Accuracy",
          pch=20,
          ylim=c(0,1), xaxt="n", outline=F,
          cex.main=2 ,
          cex.lab=2, add=T)
  

  
  
  
}


dev.off()

}










#####################################################################3
#####################################################################3
#####################################################################3

 boxplot(unlist(aff_acc_list[[d]]$random_aff_acc$RFBS2_acc),
        unlist(aff_acc_list[[d]]$random_aff_acc$DEC2_acc ), 
        unlist(aff_acc_list[[d]]$random_aff_acc$RFBS1_acc),0, 
        unlist(aff_acc_list[[d]]$random_aff_acc$RFBS0_acc), 
        unlist(aff_acc_list[[d]]$random_aff_acc$DEC0_acc ), 
        side="right", col="light gray", add=T)


{
all=c(unlist(aff_acc_list[[d]]$median_aff_acc$DEC2_acc ), 
      unlist(aff_acc_list[[d]]$median_aff_acc$RFBS2_acc))
range = c(min((all)), 
          max((all)))
plot((unlist(aff_acc_list[[d]]$median_aff_acc$DEC2_acc )) ,
     (unlist(aff_acc_list[[d]]$median_aff_acc$RFBS2_acc)),
     xlim=range, 
     ylim=range,
     xlab="log(sum(DEC true realized prob/ \n random chance true realized affinity)",
     ylab="log(sum(RFBS true realized prob/ \n random chance true realized affinity)",
     pch=20, col=tcol("black", 80)
)
abline(a = 0,b=1 )


all=c(unlist(aff_acc_list[[d]]$median_aff_acc$DEC0_acc ), 
      unlist(aff_acc_list[[d]]$median_aff_acc$RFBS0_acc))
range = c(min((all)), 
          max((all)))
plot((unlist(aff_acc_list[[d]]$median_aff_acc$DEC0_acc )) ,
     (unlist(aff_acc_list[[d]]$median_aff_acc$RFBS0_acc)),
     xlim=range, 
     ylim=range,
     xlab=" log(sum(DEC true non-affinity prob/ \n random chance true non affinity)",
     ylab="log(sum(RFBS true non-affinity prob/ \n random chance true non affinity)",
     pch=20, col=tcol("black", 80)
)

abline(a = 0,b=1 )
}


RFBS1_acc
RFBS_aff_ratios_df[[53]]
true_aff_df[[53]]

vioplot(unlist(RFBS2sim_acc), unlist(DEC2sim_acc ) )
vioplot(unlist(RFBS1sim_acc))
vioplot(unlist(RFBS0sim_acc), unlist(DEC0sim_acc ))

vioplot(unlist(RFBS2_acc), unlist(DEC2_acc )  ,side="left" )
vioplot(unlist(rRFBS2_acc), unlist(rDEC2_acc ),side="right",  add=T)

vioplot(unlist(RFBS1_acc) ,side="left"  )  
vioplot(unlist(rRFBS1_acc),side="right", add=T)

vioplot(unlist(RFBS0_acc) , unlist(DEC0_acc ) ,side="left"  )
vioplot(unlist(rRFBS0_acc), unlist(rDEC0_acc ),side="right", add=T)

vioplot(unlist(DEC2sim_acc ))
vioplot(unlist(DEC0sim_acc ))


#boxplot(unlist(RFBS_DEC0))

par(mfrow=c(2,3))

vioplot(unlist(median_RFBS2sim_acc), unlist(median_DEC2sim_acc ) )
vioplot(unlist(median_RFBS1sim_acc))
vioplot(unlist(median_RFBS0sim_acc), unlist(median_DEC0sim_acc ))

vioplot(unlist(med_RFBS2_acc),  unlist(rRFBS2_acc )  ,side="left" )
vioplot(unlist(med_DEC2_acc), unlist(rDEC2_acc ),side="right",  add=T)

vioplot(unlist(med_RFBS1_acc) ,unlist(rRFBS1_acc) )  
#vioplot(unlist(rRFBS1_acc),side="right", add=T)

vioplot(unlist(med_RFBS0_acc) , unlist(rRFBS0_acc ) ,side="left" )
vioplot(unlist(med_DEC0_acc)    , unlist(rDEC0_acc )    ,side="right", add=T)





all=c(unlist(median_DEC2sim_acc), unlist(median_RFBS2sim_acc))
range = c(min(log(all)), max(log(all)))
plot(log(unlist(median_DEC2sim_acc)) ,
     log(unlist(median_RFBS2sim_acc)),
     xlim=range, 
     ylim=range,
     xlab="log(sum(DEC true realized prob/ \n random draw of realized from statespace)",
     ylab="log(sum(RFBS true realized prob/ \n random draw of realized from any statespace)"
)
abline(a = 0,b=1 )



all=c(unlist(median_DEC0sim_acc), unlist(median_RFBS0sim_acc))
range = c(min(log(all)), max(log(all)))
plot(log(unlist(median_DEC0sim_acc)),
     log(unlist(median_RFBS0sim_acc)),
     xlim=range, ylim=range,
     xlab=" log(sum(DEC true non-affinity prob/ \n random draw of non-affinity from statespace)",
     ylab="log(sum(RFBS true non-affinity prob/ \n random draw of non-affinity from any statespace)"
)

abline(a = 0,b=1 )



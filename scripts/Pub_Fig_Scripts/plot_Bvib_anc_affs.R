{
  rm(list=ls())
  
  make_vib_clade_mtx = function(phy) {
    
    # set up clade names
    clade_mtx = matrix(NA, ncol=4, nrow=0)
    clade_mtx = rbind( clade_mtx, c( "Succotinus",    getMRCA(phy, c("V_erosum","V_betulifolium")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lobata",        getMRCA(phy, c("V_orientale","V_kansuense")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Sambucina",     getMRCA(phy, c("V_sambucinum","V_beccarii")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Coriacea",      getMRCA(phy, c("V_cylindricum","V_coriaceum")),      2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Opulus",        getMRCA(phy, c("V_opulus","V_edule")),               2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Oreinotinus",   getMRCA(phy, c("V_undulatum","V_caudatum")),         2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Dentata",       getMRCA(phy, c("V_dentatum","V_scabrellum")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Mollotinus",    getMRCA(phy, c("V_molle","V_australe")),             2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Tinus",         getMRCA(phy, c("V_tinus","V_davidii")),              2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Solenotinus",   getMRCA(phy, c("V_sieboldii","V_foetens")),          2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lutescentia",   getMRCA(phy, c("V_plicatum","V_amplifolium")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Euviburnum",    getMRCA(phy, c("V_chinshanense","V_cotinifolium")),  2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lentago",       getMRCA(phy, c("V_rufidulum","V_nudum")),            2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Punctata",      getMRCA(phy, c("V_punctatum","V_lepidotulum")),      2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Pseudotinus",   getMRCA(phy, c("V_sympodiale","V_nervosum")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Urceolata",     getMRCA(phy, c("V_urceolatum","V_taiwanianum")),     2.5, 1 ))
    clade_mtx = rbind( clade_mtx, c( "V. clemensiae", which(phy$tip.label=="V_clemensiae"),             2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "V. amplificatum", which(phy$tip.label=="V_amplificatum"),  2, 1 ))
    
    clade_mtx = rbind( clade_mtx, c( "Laminotinus",   getMRCA(phy, c("V_erosum","V_cylindricum")), 3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Porphyrotinus", getMRCA(phy, c("V_undulatum","V_australe")), 3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Crenotinus", getMRCA(phy, c("V_amplifolium","V_foetens")),   3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Valvatotinus", getMRCA(phy, c("V_lantana","V_punctatum")),   3.2, 2 ))
    
    colnames(clade_mtx) = c("name","node","level","barsize")
    clade_df = data.frame(clade_mtx, stringsAsFactors=F)
    clade_df$node = as.numeric(clade_df$node)
    clade_df$level = as.numeric(clade_df$level)
    clade_df$barsize = as.numeric(clade_df$barsize)
    
    return(clade_df)
    
  }
  
  circ_node_coords=function (tr, coords_fun = "plot_phylo3_nodecoords", tmplocation = "auto", 
                             root.edge = TRUE) 
  {
    defaults = "\n\tcoords_fun=\"plot_phylo3_nodecoords\"\n\ttmplocation=\"manual\"\n\t"
    if (tmplocation == "manual") {
      scriptdir = "/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/"
    }
    else if (tmplocation == "auto") {
      extdata_dir = np(system.file("extdata", package = "BioGeoBEARS"))
      scriptdir = slashslash(paste(extdata_dir, "a_scripts", 
                                   sep = "/"))
    }
    else {
      scriptdir = tmplocation
    }
    file_to_source = np(paste(addslash(scriptdir), coords_fun, 
                              ".R", sep = ""))
    source(file = file_to_source)
    trcoords = matrix(data = c(1, 1), nrow = 1, ncol = 2)
    if (packageVersion("ape") < 5) {
      cmdstr = paste("trcoords = ", coords_fun, "(tr, plot=FALSE, root.edge=root.edge)", 
                     sep = "")
      eval(parse(text = cmdstr))
      nodecoords = cbind(trcoords$xx, trcoords$yy)
    }
    else {
      trcoords = plot_phylo3_nodecoords_APE5(tr, plot = FALSE, 
                                             root.edge = root.edge,type = "fan")
      nodecoords = cbind(trcoords$xx, trcoords$yy)
    }
    nodes_xy = adf(nodecoords)
    names(nodes_xy) = c("x", "y")
    row.names(nodes_xy) = 1:nrow(nodecoords)
    return(nodes_xy)
  }
  
  
  plotRFBSancaff=function(tree, state_mat, aff_ratios, plot_name=NULL, clade_labels=NULL, dims=c(80,150), plot=c() ){
    
    fan=F
    nbiomes=ncol(aff_ratios$node_fund)
    nstates=nrow(state_mat)
    node_length=nrow(aff_ratios$node_non)
    
    no_plot_children=no_plot-length(tree$tip.label)
    plot_children=(plot[-(1:length(tree$tip.label))]-length(tree$tip.label))
    
    sq = seq(-nbiomes*0.3,nbiomes*0.3, length.out = nbiomes)
    
    
    
    cexcircles = (sq[2]- sq[1])*.5
    
    #sq=c(0,0,0)
    
    tip_offset=sum(abs(range(sq )))*2
    
    
    
    #####get colors#####
    library(RColorBrewer)
    library(BioGeoBEARS)
    library(yarrr)
    n <- nstates
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector_full= unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector=col_vector_full[1:n]
    
    if(fan==T){
      
      coords=circ_corner_coords(tree)
      
    }else{
      
      coords=BioGeoBEARS::corner_coords(tree,tmplocation="auto")
    }
    
    
    #####get pie charts for nodes##### 
    
    # sort probabilities for corners
    #find right children nodes
    right_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
    #right_child_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
    
    right_child_fund_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
    right_child_real_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
    
    #sim_right_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
    
    #tree$edge prob, from what i can tell the topmost node connection in the matrix is for the right child and the second for the left
    for (i in 1:length(unique(tree$edge[,1]))){
      
      nd=i+length(tree$tip.label)  
      row_index <- which(tree$edge[,1]==nd, arr.ind = TRUE)[1]
      child <- tree$edge[row_index,2]
      
      #right_child_probs[i,] = corner_state_freqs[child,]
      #sim_right_child_probs[i,] = sim_corner_state_freqs[child,]
      #right_child_probs[i,] = corner_state_freqs[child,]
      right_child_fund_af_ratio[i,] = aff_ratios$corner_fund[child,]
      right_child_real_af_ratio[i,] = aff_ratios$corner_real[child,]
      
    }
    
    aff_ratios
    #find left children nodes
    
    left_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
    #sim_left_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
    left_child_fund_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
    left_child_real_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
    
    
    for (i in 1:length(unique(tree$edge[,1]))){
      
      nd=i+length(tree$tip.label)  
      row_index <- which(tree$edge[,1]==nd, arr.ind = TRUE)[2]
      child <- tree$edge[row_index,2]
      
     # left_child_probs[i,] = corner_state_freqs[child,]
      #sim_left_child_probs[i,] = sim_corner_state_freqs[child,]
      
      left_child_fund_af_ratio[i,] = aff_ratios$corner_fund[child,]
      left_child_real_af_ratio[i,] = aff_ratios$corner_real[child,]
      
    }
    #points(coords$leftcorns[,2:3])
    
    ######plot #####  
    
    
    #st2leg=colSums(node_state_freqs)>1.0
    
    
    
    
    
    
    
    ###### ############
    
    
    
    #node_chords=circ_node_coords(tree)
    
    if(fan==T){
      
      node_chords=circ_node_coords(tree)
      
    }else{
      
      node_chords=node_coords(tree)
    }
    
    
    {
      
      if(!is.null(plot_name)){
        pdf(plot_name,width = dims[[1]],height = dims[[2]])  
      
      }
      plot(tree,label.offset=tip_offset, cex=3) 
      #plot_phylo3_nodecoords_APE5(tree, plot = T, 
      #                            root.edge = F, type = "fan")
      #
      for (nd in (1:node_length)[plot]){
        
        if(nd>length(tree$tip.label)){
          
          
          floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red")  ,  radius =cexcircles, border=c("white", "red")  )
          floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles, border=c("white", "green"))
          floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue") , radius  =cexcircles, border=c("white", "blue") )
          
          floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "red")  )
          floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "green"))
          floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "blue") )
          
          #points(node_chords[nd,1]+sq[[1]], node_chords[nd,2], pch=19,col=transparent("red"  , trans.val = 1-aff_ratios$node_fund[nd,1]), cex=cexcircles)
          #points(node_chords[nd,1]+sq[[2]], node_chords[nd,2], pch=19,col=transparent("green"  , trans.val = 1-aff_ratios$node_fund[nd,2]), cex=cexcircles)
          #points(node_chords[nd,1]+sq[[3]], node_chords[nd,2], pch=19,col=transparent("light blue" , trans.val = 1-aff_ratios$node_fund[nd,3]), cex=cexcircles)
          #
          #points(node_chords[,1]+sq[[1]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          #points(node_chords[,1]+sq[[2]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          #points(node_chords[,1]+sq[[3]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          
          floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5,border=c("white", "red")  )
          floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5,border=c("white", "green"))
          floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5,border=c("white", "blue") )
          
          #points(node_chords[nd,1]+sq[[1]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,1]), cex=cexcircles*.5)
          #points(node_chords[nd,1]+sq[[2]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,2]), cex=cexcircles*.5)
          #points(node_chords[nd,1]+sq[[3]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,3]), cex=cexcircles*.5)
          
        }else{
          #points(node_chords[nd,1]+sq[[1]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("red"  , trans.val = 1-aff_ratios$node_fund[nd,1]), cex=cexcircles)
          #points(node_chords[nd,1]+sq[[2]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("green", trans.val = 1-aff_ratios$node_fund[nd,2]), cex=cexcircles)
          #points(node_chords[nd,1]+sq[[3]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("light blue" , trans.val = 1-aff_ratios$node_fund[nd,3]), cex=cexcircles)
          #  points(node_chords[,1]+sq[[1]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          #  points(node_chords[,1]+sq[[2]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          #  points(node_chords[,1]+sq[[3]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
          #points(node_chords[nd,1]+sq[[1]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,1]), cex=cexcircles*.5)
          #points(node_chords[nd,1]+sq[[2]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,2]), cex=cexcircles*.5)
          #points(node_chords[nd,1]+sq[[3]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,3]), cex=cexcircles*.5)
          
          
          if(nd%%1==0){
            
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles  , border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles, border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles , border=c("white", "blue") )
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "blue") )
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "blue") )
          }else{
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles  , border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles, border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles , border=c("white", "blue") )
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5,  border=c("white", "blue") )
            floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "red")  )
            floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "green"))
            floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5, border=c("white", "blue") )
            
            
            
            
            
          }
          
        }
        
      }
      
      for(nd in (1:nrow(coords$rightcorns))[plot_children]){
        
        
        
        
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,1],  right_child_fund_af_ratio[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles   ,border=c("white", "red")  )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,2],  right_child_fund_af_ratio[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles ,border=c("white", "green"))
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,3],  right_child_fund_af_ratio[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles  ,border=c("white", "blue") )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5, border=c("white", "red")  )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5, border=c("white", "green"))
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5, border=c("white", "blue") )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,1],  right_child_real_af_ratio[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*.5,border=c("white", "red")   )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,2],  right_child_real_af_ratio[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*.5,border=c("white", "green") )
        floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,3],  right_child_real_af_ratio[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*.5,border=c("white", "blue")  )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,1],  left_child_fund_af_ratio[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles   ,border=c("white", "red")  )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,2],  left_child_fund_af_ratio[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles ,border=c("white", "green"))
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,3],  left_child_fund_af_ratio[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles  ,border=c("white", "blue") )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5 ,border=c("white", "red")  )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5 ,border=c("white", "green"))
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5 ,border=c("white", "blue") )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,1],  left_child_real_af_ratio[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*.5, border=c("white", "red")  )
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,2],  left_child_real_af_ratio[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*.5, border=c("white", "green"))
        floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,3],  left_child_real_af_ratio[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*.5, border=c("white", "blue") )
        
        
        #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-right_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
        #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-  right_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
        #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1- right_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
        #
        #
        #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("red", trans.val = 1- left_child_fund_af_ratio[nd,1]), cex=cexcircles)
        #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("green"  , trans.val = 1-   left_child_fund_af_ratio[nd,2]), cex=cexcircles)
        #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("light blue" , trans.val = 1-  left_child_fund_af_ratio[nd,3]), cex=cexcircles)
        #  points(coords$leftcorns[,2]+sq[[1]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #  points(coords$leftcorns[,2]+sq[[2]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #  points(coords$leftcorns[,2]+sq[[3]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-left_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
        #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-  left_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
        #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1- left_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
        
        #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("red", trans.val = 1- right_child_fund_af_ratio[nd,1]), cex=cexcircles)
        #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("green"  , trans.val =   1- right_child_fund_af_ratio[nd,2]), cex=cexcircles)
        #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("light blue"  , trans.val =  1- right_child_fund_af_ratio[nd,3]), cex=cexcircles)
        #points(coords$rightcorns[,2]+sq[[1]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$rightcorns[,2]+sq[[2]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$rightcorns[,2]+sq[[3]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-right_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
        #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-  right_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
        #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1- right_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
        #
        #
        #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("red", trans.val = 1- left_child_fund_af_ratio[nd,1]), cex=cexcircles)
        #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("green"  , trans.val = 1-   left_child_fund_af_ratio[nd,2]), cex=cexcircles)
        #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("light blue" , trans.val = 1-  left_child_fund_af_ratio[nd,3]), cex=cexcircles)
        #points(coords$leftcorns[,2]+sq[[1]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$leftcorns[,2]+sq[[2]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$leftcorns[,2]+sq[[3]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
        #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-left_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
        #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-  left_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
        #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1- left_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
        
        
      }
      
      #nodelabels()
      
      
      if(!is.null(clade_labels)){
      for ( i in 1:nrow(clade_labels)){
        
        cladelabels(tree =tree
                    ,node   = clade_labels$node[[i]]
                    ,offset =clade_labels$level[[i]]*7
                    , text  = clade_labels$name[[i]] 
                    ,wing.length = 0,
                     cex=3 )
      }
      }
      
      
      if(!is.null(plot_name)){
        dev.off()        
      }
      
      
    }
    
    
  }
  
  
  emp_anc_affs_from_dir=function(dir_name, states_mat){
    
    file_list=list.files(dir_name)
    run_list=file_list[grep("_log_anc_states",unlist(file_list),fixed=FALSE)]
    clado_list=file_list[grep("_log_anc_clados",unlist(file_list),fixed=FALSE)]
    
    
    nbiomes=ncol(states_mat)
    #states_list=unlist(lapply(1:nrow(states_mat), function(i) paste(states_mat[i,],collapse = " ")))
    anc_states_string_vec=list()
    anc_states=list()
    burnin=.5
    
    ####get anc state chains and modes #######
    {
      
      anc_states_chains=do.call(rbind, lapply(run_list, function(run) read.table(paste(dir_name,run, sep="/"), header = F)))
      
      anc_clado_chains=do.call(rbind, lapply(clado_list, function(run) read.table(paste(dir_name,run, sep="/"), header = F)))
      
      node_length=ncol(anc_states_chains)/2
      
      burnin=nrow(anc_states_chains)/2
      
    }
    
    
    anc_states_post_burn= do.call(cbind, lapply(1:ncol(anc_states_chains), function(i)  anc_states_chains[ burnin:nrow( anc_states_chains),i] ))
    
    aff_ratios=get_aff_ratios(anc_states_post_burn, states_mat)
    
    return(aff_ratios)
    
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

setwd("/Volumes/michael.landis/Active/RFBS_RIS")

workdir=getwd()


dir_names=c("Bvib_3nB_LN_0_p5_1000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_LN_0_p5_1000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_LN_0_p5_1000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_LN_0_p5_1000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_LN_0_p5_1000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_LN_0_p5_1000000_DEC_eco_clado_2g_2l_gl")
            
dir_names=c(
            "Bvib_3nB_Exp0p5_1000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_DEC_eco_clado_2g_2l_gl",
            "Bvib_3nB_Exp0p5_1000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_admat_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
            "Bvib_3nB_Exp0p5_1000000_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds"
           
)

new_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_3000000",list.files()  ))]
#old_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_10000000",list.files()  ))][c( 5,1,6,7, 2,3,8,4)]


dir_names=c(old_RFBS_dir_names, new_RFBS_dir_names)

dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_3000000",list.files()  )[grep("Bvib_3nB_Exp0p5_3000000",list.files()  ) %in% grep("Foss",list.files()  )])]


burnin=0.5
RFBS_states_space= read.table("rf_states.txt")
DEC_states_space=read.table("DEC_states.txt")

tree=read.tree("viburnum_data_files/viburnum_sorted.tre")


states_mat=read.table(paste("rf_states.txt", sep="/"))
#states_mat=read.table(paste("DEC_states.txt", sep="/"))



affs_list=list()
length(dir_names)
dir_inds=c(14,28, 29, 42, 57, 71, 85, 99, 113 )

for(d in dir_inds){
  print(d)
  if(grepl("DEC",dir_names[[d]])){
    
    affs_list[[d]]=emp_anc_affs_from_dir(dir_names[[d]], DEC_states_space)
    
    
  }
  
  affs_list[[d]]=emp_anc_affs_from_dir(dir_names[[d]], RFBS_states_space)
}


affs_list[[1]]

no_plot=c(175:215, 217:219, 223:225, 226:239)
no_plot_children=no_plot-length(tree$tip.label)

plot=c(1:168, 169:173, 174, 202, 204, 208, 209, 210, 
       211, 212, 216:219, 221, 231, 244, 250, 251:252, 256, 257,  
       261,265, 264, 266,  267, 272:277, 279, 284,286, 290, 292:295,297:299, 302,304, 306, 309, 311, 314,315, 321:324)
plot_children=(plot[-(1:168)]-length(tree$tip.label))

plot_name = paste(dir_names[[d]],"anc_aff.pdf", sep="/")

plot_names=NULL
dir_names

clade_labels=make_vib_clade_mtx(tree)


{
  pdf("tree.pdf", width = 100, height = 100)
  
  plot(tree)
  
  nodelabels(cex=5)
  
  dev.off()
}


for (d in dir_inds){
  print(d)
    
  plot_name=paste("rf_sim_plots/", dir_names[[d]],"_anc_aff.pdf", sep="")
  
  plotRFBSancaff(tree, states_mat,affs_list[[d]], plot_name,clade_labels, dims=c(80,150), plot)
}



########################################################################
plot(tree)
for ( i in 1:nrow(clade_labels)){

cladelabels(tree =tree
            ,node   = clade_labels$node[[i]]
            ,offset =clade_labels$level[[i]]
            , text  = clade_labels$name[[i]] ,
            ,wing.length = 0)
}
clade_labels=make_vib_clade_mtx(tree)


segments(tree)



pp$edge

points(pp$xx[[162]], pp$yy[[162]], col=2)
points(0,0,col=2)
pp$edge



coord_rotate = function(x, angle=0, origin=c(0,0)) {
 
    xy = x
    xy_new = xy
    xy_new[1] = origin[1] + cos(angle) * (xy[1] - origin[1]) - sin(angle) * (xy[2] - origin[2])
    xy_new[2] = origin[2] + sin(angle) * (xy[1] - origin[1]) + cos(angle) * (xy[2] - origin[2])
    x = xy_new

  return(x)
}







plot(tree,plot=FALSE, type="fan") 

pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

get_phylo_fan_coords=function(tree, pp=NULL){
  
  if(is.null(pp)==T){
    
    plot(tree,plot=FALSE, type="fan") 
    
    pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    
  }

node_coords=cbind(1:length(pp$xx), pp$xx, pp$yy)
rightcorn_coords=node_coords
leftcorn_coords=node_coords

for (i in unique(pp$edge[,1])){

  node=i
  
  children=pp$edge[pp$edge[,1]==node,2]
  
  for (x in c(1,2)){
  
    c=children[[x]]
    P1=c(0,0)
    P2=c(pp$xx[[node]], pp$yy[[node]])
    P3=c(pp$xx[[c]], pp$yy[[c]])
    
    
    result = atan2(P3[[2]] - P1[[2]], P3[[1]] - P1[[1]]) -
      atan2(P2[[2]] - P1[[2]], P2[[1]] - P1[[1]])
    
    rotated_pt=coord_rotate(P2, angle=result)
    
    if(x==1){
     
      rightcorn_coords[rightcorn_coords[,1]==i,2:3] =  rotated_pt
      
    }else{
      
      leftcorn_coords[leftcorn_coords[,1]==i,2:3]  =  rotated_pt
    
    }
    #points(x = rotated_pt[[1]], y=rotated_pt[[2]], col=3)
  }
}

return(list("nodes"=node_coords,
            "right"=rightcorn_coords,
            "left"= leftcorn_coords))

}

get_vec_angle=function(P1,P2,P3){
  #P1 is origin of both vectors
  if("numeric" %in% class(P3)){
  
  result = atan2(P3[[2]] - P1[[2]], P3[[1]] - P1[[1]]) -
    atan2(P2[[2]] - P1[[2]], P2[[1]] - P1[[1]])
  }else{
    result = atan2(P3[,2] - P1[[2]], P3[,1] - P1[[1]]) -
      atan2(P2[,2] - P1[[2]], P2[,1] - P1[[1]])    
  }
  
  return(result)
}

move_coord_by_angle=function(x,y, ang_deg, distance){
 
  
  ang_rad = ang_deg #* pi / 180
  x1 = x + sin(ang_rad)* distance
y1 = y + cos(ang_rad) * distance

return(cbind(x1, y1))
}


move_coord_from_origin=function(coord, origin){
  
  
  
  
  
  
}
  




plot(tree, type="fan") 
points(rightcorn_coords[,2:3], col=2)
points(leftcorn_coords[,2:3], col=2)
points(node_coords[,2:3], col=2)

{
tree_height=max(nodeHeights(tree))
pie_dist_mod=0.05
biome_aff_dist=tree_height*pie_dist_mod
nbiomes=3
pie_spaces= seq(-nbiomes*biome_aff_dist,nbiomes*biome_aff_dist, length.out = nbiomes)

scalar=1

phy_coords=get_phylo_fan_coords(tree)


origin_angles=atan2(phy_coords$nodes[,3], phy_coords$nodes[,2])

origin_angles[origin_angles<0]=origin_angles[origin_angles<0]+2*pi

#origin_angles=get_vec_angle(c(0,0),cbind(phy_coords$nodes[,2],0),phy_coords$nodes[,2:3])

phy_coords$left[,2:3]



phy_coords$nodes[,2]/phy_coords$nodes[,3]


coord_rotate()

rotate_pt_around_origin=function(origin, start_angle_pt, end_angle_pt){

P1=c(0,0)
P2=c(pp$xx[[node]], pp$yy[[node]])
P3=c(pp$xx[[c]], pp$yy[[c]])


result = atan2(P3[[2]] - P1[[2]], P3[[1]] - P1[[1]]) -
  atan2(P2[[2]] - P1[[2]], P2[[1]] - P1[[1]])

rotated_pt=coord_rotate(P2, angle=result)

}

biome_node_coords=lapply(pie_spaces, function(b) move_coord_by_angle(phy_coords$nodes[,2],
                                                                     phy_coords$nodes[,3], origin_angles[[1]], b))

biome_leftcorn_coords=lapply(pie_spaces, function(b) move_coord_by_angle(phy_coords$left[,2], 
                                                                         phy_coords$left[,3], 
                                                                         origin_angles, b))
biome_rightcorn_coords=lapply(pie_spaces, function(b) move_coord_by_angle(phy_coords$right[,2], 
                                                                          phy_coords$right[,3], 
                                                                          origin_angles, b))


plot(tree, type="fan")
length(pie_spaces)
for(b in 2:3){
points(biome_node_coords[[b]], col=b)
}
  

}



node_coords[,2:3]

x=node_coords[1,2]
y=node_coords[1,2]


if(x<0){
  
  
}


points(leftcorn_coords[,2:3], col=2)
points(node_coords[,2:3], col=2)



###############################################

plotRFBSancaff_fan=function(tree, state_mat, aff_ratios, plot_name=NULL, clade_labels=NULL, dims=c(80,150), plot=c() ){
  
  fan=F
  nbiomes=ncol(aff_ratios$node_fund)
  nstates=nrow(state_mat)
  node_length=nrow(aff_ratios$node_non)
  
  no_plot_children=no_plot-length(tree$tip.label)
  plot_children=(plot[-(1:161)]-length(tree$tip.label))
  
  sq = seq(-nbiomes*0.3,nbiomes*0.3, length.out = nbiomes)
  
  
  
  cexcircles = (sq[2]- sq[1])*.5
  
  #sq=c(0,0,0)
  
  tip_offset=sum(abs(range(sq )))*2
  
  
  
  #####get colors#####
  library(RColorBrewer)
  library(BioGeoBEARS)
  library(yarrr)
  n <- nstates
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector_full= unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector=col_vector_full[1:n]
  
  if(fan==T){
    
    coords=circ_corner_coords(tree)
    
  }else{
    
    coords=BioGeoBEARS::corner_coords(tree,tmplocation="auto")
  }
  
  
  #####get pie charts for nodes##### 
  
  # sort probabilities for corners
  #find right children nodes
  right_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
  #right_child_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
  
  right_child_fund_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
  right_child_real_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
  
  #sim_right_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
  
  #tree$edge prob, from what i can tell the topmost node connection in the matrix is for the right child and the second for the left
  for (i in 1:length(unique(tree$edge[,1]))){
    
    nd=i+length(tree$tip.label)  
    row_index <- which(tree$edge[,1]==nd, arr.ind = TRUE)[1]
    child <- tree$edge[row_index,2]
    
    #right_child_probs[i,] = corner_state_freqs[child,]
    #sim_right_child_probs[i,] = sim_corner_state_freqs[child,]
    #right_child_probs[i,] = corner_state_freqs[child,]
    right_child_fund_af_ratio[i,] = aff_ratios$corner_fund[child,]
    right_child_real_af_ratio[i,] = aff_ratios$corner_real[child,]
    
  }
  
  aff_ratios
  #find left children nodes
  
  left_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
  #sim_left_child_probs=matrix(0, ncol=nstates,nrow=(node_length-1)/2)
  left_child_fund_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
  left_child_real_af_ratio=matrix(0, ncol=nbiomes,nrow=(node_length-1)/2)
  
  
  for (i in 1:length(unique(tree$edge[,1]))){
    
    nd=i+length(tree$tip.label)  
    row_index <- which(tree$edge[,1]==nd, arr.ind = TRUE)[2]
    child <- tree$edge[row_index,2]
    
    # left_child_probs[i,] = corner_state_freqs[child,]
    #sim_left_child_probs[i,] = sim_corner_state_freqs[child,]
    
    left_child_fund_af_ratio[i,] = aff_ratios$corner_fund[child,]
    left_child_real_af_ratio[i,] = aff_ratios$corner_real[child,]
    
  }
  #points(coords$leftcorns[,2:3])
  
  ######plot #####  
  
  
  #st2leg=colSums(node_state_freqs)>1.0
  
  
  
  
  
  
  
  ###### ############
  
  
  
  #node_chords=circ_node_coords(tree)
  
  if(fan==T){
    
    node_chords=circ_node_coords(tree)
    
  }else{
    
    node_chords=node_coords(tree)
  }
  
  
  {
    
    if(!is.null(plot_name)){
      pdf(plot_name,width = dims[[1]],height = dims[[2]])  
      
    }
    plot(tree,label.offset=tip_offset, cex=3) 
    #plot_phylo3_nodecoords_APE5(tree, plot = T, 
    #                            root.edge = F, type = "fan")
    #
    for (nd in (1:node_length)[plot]){
      
      if(nd>length(tree$tip.label)){
        
        
        floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles)
        floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles)
        floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles)
        
        floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
        floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
        floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
        
        #points(node_chords[nd,1]+sq[[1]], node_chords[nd,2], pch=19,col=transparent("red"  , trans.val = 1-aff_ratios$node_fund[nd,1]), cex=cexcircles)
        #points(node_chords[nd,1]+sq[[2]], node_chords[nd,2], pch=19,col=transparent("green"  , trans.val = 1-aff_ratios$node_fund[nd,2]), cex=cexcircles)
        #points(node_chords[nd,1]+sq[[3]], node_chords[nd,2], pch=19,col=transparent("light blue" , trans.val = 1-aff_ratios$node_fund[nd,3]), cex=cexcircles)
        #
        #points(node_chords[,1]+sq[[1]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        #points(node_chords[,1]+sq[[2]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        #points(node_chords[,1]+sq[[3]], node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        
        floating.pie.asp(node_chords[nd,1]+sq[[1]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
        floating.pie.asp(node_chords[nd,1]+sq[[2]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
        floating.pie.asp(node_chords[nd,1]+sq[[3]], node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
        
        #points(node_chords[nd,1]+sq[[1]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,1]), cex=cexcircles*.5)
        #points(node_chords[nd,1]+sq[[2]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,2]), cex=cexcircles*.5)
        #points(node_chords[nd,1]+sq[[3]], node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,3]), cex=cexcircles*.5)
        
      }else{
        #points(node_chords[nd,1]+sq[[1]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("red"  , trans.val = 1-aff_ratios$node_fund[nd,1]), cex=cexcircles)
        #points(node_chords[nd,1]+sq[[2]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("green", trans.val = 1-aff_ratios$node_fund[nd,2]), cex=cexcircles)
        #points(node_chords[nd,1]+sq[[3]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent("light blue" , trans.val = 1-aff_ratios$node_fund[nd,3]), cex=cexcircles)
        #  points(node_chords[,1]+sq[[1]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        #  points(node_chords[,1]+sq[[2]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        #  points(node_chords[,1]+sq[[3]]+tip_offset/2, node_chords[,2], pch=19,col=0, cex=cexcircles*.6)
        #points(node_chords[nd,1]+sq[[1]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,1]), cex=cexcircles*.5)
        #points(node_chords[nd,1]+sq[[2]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,2]), cex=cexcircles*.5)
        #points(node_chords[nd,1]+sq[[3]]+tip_offset/2, node_chords[nd,2], pch=19,col=transparent(1, trans.val = 1-aff_ratios$node_real[nd,3]), cex=cexcircles*.5)
        
        
        if(nd%%1==0){
          
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(1/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
        }else{
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,1],  aff_ratios$node_fund[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,2],  aff_ratios$node_fund[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1- aff_ratios$node_fund[nd,3],  aff_ratios$node_fund[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles)
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1), pch=19,col=c(0), radius  =cexcircles*.5)
          floating.pie.asp(node_chords[nd,1]+sq[[1]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,1],   aff_ratios$node_real[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
          floating.pie.asp(node_chords[nd,1]+sq[[2]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,2],   aff_ratios$node_real[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
          floating.pie.asp(node_chords[nd,1]+sq[[3]]+tip_offset*(2/3), node_chords[nd,2], x=c(1-aff_ratios$node_real[nd,3],   aff_ratios$node_real[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*0.5)
          
          
          
          
          
        }
        
      }
      
    }
    
    for(nd in (1:nrow(coords$rightcorns))[plot_children]){
      
      
      
      
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,1],  right_child_fund_af_ratio[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,2],  right_child_fund_af_ratio[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_fund_af_ratio[nd,3],  right_child_fund_af_ratio[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,1],  right_child_real_af_ratio[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,2],  right_child_real_af_ratio[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      floating.pie.asp(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.4, x=c(1- right_child_real_af_ratio[nd,3],  right_child_real_af_ratio[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,1],  left_child_fund_af_ratio[nd,1]), pch=19,col=c("white", "red"), radius  =cexcircles)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,2],  left_child_fund_af_ratio[nd,2]), pch=19,col=c("white", "green"), radius  =cexcircles)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_fund_af_ratio[nd,3],  left_child_fund_af_ratio[nd,3]), pch=19,col=c("white", "blue"), radius  =cexcircles)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1), pch=19,col=0, radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,1],  left_child_real_af_ratio[nd,1]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,2],  left_child_real_af_ratio[nd,2]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      floating.pie.asp(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.4, x=c(1- left_child_real_af_ratio[nd,3],  left_child_real_af_ratio[nd,3]), pch=19,col=c("white", "black"), radius  =cexcircles*.5)
      
      
      #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-right_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
      #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-  right_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
      #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1- right_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
      #
      #
      #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("red", trans.val = 1- left_child_fund_af_ratio[nd,1]), cex=cexcircles)
      #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("green"  , trans.val = 1-   left_child_fund_af_ratio[nd,2]), cex=cexcircles)
      #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("light blue" , trans.val = 1-  left_child_fund_af_ratio[nd,3]), cex=cexcircles)
      #  points(coords$leftcorns[,2]+sq[[1]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #  points(coords$leftcorns[,2]+sq[[2]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #  points(coords$leftcorns[,2]+sq[[3]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-left_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
      #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-  left_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
      #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1- left_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
      
      #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("red", trans.val = 1- right_child_fund_af_ratio[nd,1]), cex=cexcircles)
      #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("green"  , trans.val =   1- right_child_fund_af_ratio[nd,2]), cex=cexcircles)
      #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent("light blue"  , trans.val =  1- right_child_fund_af_ratio[nd,3]), cex=cexcircles)
      #points(coords$rightcorns[,2]+sq[[1]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$rightcorns[,2]+sq[[2]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$rightcorns[,2]+sq[[3]],   coords$rightcorns[,3]-0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$rightcorns[nd,2]+sq[[1]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-right_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
      #points(coords$rightcorns[nd,2]+sq[[2]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1-  right_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
      #points(coords$rightcorns[nd,2]+sq[[3]], coords$rightcorns[nd,3]-0.1, pch=19,col=transparent(1, trans.val = 1- right_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
      #
      #
      #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("red", trans.val = 1- left_child_fund_af_ratio[nd,1]), cex=cexcircles)
      #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("green"  , trans.val = 1-   left_child_fund_af_ratio[nd,2]), cex=cexcircles)
      #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent("light blue" , trans.val = 1-  left_child_fund_af_ratio[nd,3]), cex=cexcircles)
      #points(coords$leftcorns[,2]+sq[[1]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$leftcorns[,2]+sq[[2]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$leftcorns[,2]+sq[[3]],   coords$leftcorns[,3]+0.1, pch=19,col=0, cex=cexcircles*.6)
      #points(coords$leftcorns[nd,2]+sq[[1]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-left_child_real_af_ratio[nd,1]), cex=cexcircles*.5)
      #points(coords$leftcorns[nd,2]+sq[[2]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1-  left_child_real_af_ratio[nd,2]), cex=cexcircles*.5)
      #points(coords$leftcorns[nd,2]+sq[[3]], coords$leftcorns[nd,3]+0.1, pch=19,col=transparent(1, trans.val = 1- left_child_real_af_ratio[nd,3]), cex=cexcircles*.5)
      
      
    }
    
    #nodelabels()
    
    
    if(!is.null(clade_labels)){
      for ( i in 1:nrow(clade_labels)){
        
        cladelabels(tree =tree
                    ,node   = clade_labels$node[[i]]
                    ,offset =clade_labels$level[[i]]*7
                    , text  = clade_labels$name[[i]] 
                    ,wing.length = 0,
                    cex=3 )
      }
    }
    
    
    if(!is.null(plot_name)){
      dev.off()        
    }
    
    
  }
  
  
}



library(DELTD)

# Number of replicates
n = 99

# HPD coverage
p = 0.5

#calc_HPD_ci(100, 0.95)

calc_HPD_ci=function(n, p){
  # Mean, E[X]
  mu = n * p
  
  # Var, E[X^2]
  variance = n * p * (1-p)
  
  # stddev(X)
  sigma = sqrt(variance)
  
  # 95% of similar experiments would produce coverage frequencies
  # within the confidence interval, mu +/- 1.96 * sigma
  ci = (mu + 1.96 * c(-sigma, +sigma))/n
  
  return(list("ci"=ci,
              "sigma"=sigma,
              "var"=variance,
              "mean"=mu))
  
}

calc_HPD_ci(100, 0.95)


{
  library(stringr)
  library(dplyr)
  library(coda)
  library(phytools)
  
  
  createLayoutMatrix <- function(Nrow, Ncol, blockRows, blockCols, blockFillOrder = "byrow", gridFillOrder = "byrow") {
    # Create an empty matrix for the layout
    layoutMatrix <- matrix(0, nrow = Nrow * blockRows, ncol = Ncol * blockCols)
    
    # Initialize the plot number
    plotNumber <- 1
    
    # Define a function to increment the plot number
    incrementPlotNumber <- function() {
      num <- plotNumber
      plotNumber <<- plotNumber + 1
      return(num)
    }
    
    # Fill the layout matrix
    for (gridIndex in 1:(blockRows * blockCols)) {
      # Determine the block's position in the grid
      if (gridFillOrder == "byrow") {
        gridRow <- (gridIndex - 1) %/% blockCols
        gridCol <- (gridIndex - 1) %% blockCols
      } else {
        gridRow <- (gridIndex - 1) %% blockRows
        gridCol <- (gridIndex - 1) %/% blockRows
      }
      
      # Fill the block
      for (blockIndex in 1:(Nrow * Ncol)) {
        if (blockFillOrder == "byrow") {
          blockRow <- (blockIndex - 1) %/% Ncol
          blockCol <- (blockIndex - 1) %% Ncol
        } else {
          blockRow <- (blockIndex - 1) %% Nrow
          blockCol <- (blockIndex - 1) %/% Nrow
        }
        
        # Calculate the position in the layout matrix
        row <- gridRow * Nrow + blockRow + 1
        col <- gridCol * Ncol + blockCol + 1
        
        # Place the plot number
        layoutMatrix[row, col] <- incrementPlotNumber()
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

    dir_names=
      c(        "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_1_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_rf_gl_ds",
                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_1_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_1_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_1_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
                "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_1_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds")
    
    dir_names=
    c("BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds"
    )
    
    
   #dir_names = 
   #  c(  "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_rf_gl_ds" ,
   #     "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds",
   #     "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds",
   #     "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds",
   #     "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds"
   #     )
        
        dir_names=
    c("BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_2sw_dr_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.25]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.33][0.75]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.25]_rf_gl_ds"
   #,"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_5000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_1g_1l_dr_unce[1.0][0.66][0.75]_rf_gl_ds"
    )


    dir_labels=c("Full data",
                 #"50 Tips / 1.0 rate prior",
                 "0.33 biomes / 0.25 tips",
                 "0.33 biomes / 0.75 tips",
                 #"150 Tips / 1.0 rate prior",
                 #"150 Tips / 5.0 rate prior",
                 "0.66 biomes / 0.25 tips",
                 "0.66 biomes / 0.75 tips"
)
    
    
    dir_labels=c("Full data"
    )
    
    #dir_names=c("BSim_RFBS_DEC_compAbs150t_3nB_LN_0_0p5iter_5000_clado_Rtre_2g_1g_1l_dr")
    #dir_names=c("saved_BSim_runs/BSim_RFBS_DEC_compAbs150t_3nB_LN_0_0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
    #                     "[1.0][0.25][0.66]",
    #"[1.0][0.5][0.33]",
    #"[1.0][0.5][0.66]",
    #"[1.0][1.0][1.0]")
    
    #dir_names=c(#"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",   
    #            #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
    
    
    #dir_names=c(#"BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_1g_1l_dr_rf_gl",
    #           # "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_rf_ds",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p1iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
    #            #"BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_1g_1l_dr_rf_gl",
    #           # "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf_ds",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp0p5iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds",
    #            #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_1g_1l_dr_rf_gl",
    #            #"BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_gl",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_ds",
    #            "BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_300000_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")
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
    #     
    #          
    #   # dir_names=c("saved_BSim_runs/Draft_2/BSim_RFBS_DEC_compPO_50t_3nB_Exp0p5iter_500000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds")           
    #            
    #dir_labels=c("no missing data",
    #             "-2 ambiguous states in 25% tips",
    #             "-2 ambiguous states in 50% tips",
    #             "-2 ambiguous states in 75% tips",
    #             "-2 ambiguous states in 100% tips",
    #             "-1 ambiguous states in 25% tips",
    #             "-1 ambiguous states in 50% tips",
    #             "-1 ambiguous states in 75% tips",
    #             "-1 ambiguous states in 100% tips",
    #             "all ambiguous states in 100% tips")
    #
    # 
    
    # dir_ind=c(1:10)
    
    
    post_median_plots <- vector(length(dir_names), mode='list')
    
    
    plot_prior=T
    
    prior_dir="BSim_RFBS_DEC_compPO_50t_3nB_Exp0p5iter_500000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds"
    #if running files from HPC
    setwd("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim")
    
    workdir=getwd()
    
    #dir_name=dir_names[[1]]
    
    sim_dat=T
    
    library(stringr)
    
    library(dplyr)
    library(coda)
    
    
    
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
    
  }
  






sim_pars_total_list=list()
post_median_total_list=list()
lower_HPD_list=list()
upper_HPD_list=list()
isin_HPD_list=list()
chain_list_list=list()
test_list=list()
isinHPD=list()
switch_on_list = list()
switch_on_PD_list = list()
correct_switch_ratio_list = list()
savage_dickey_ratios_list=list()


log_dirs     =  lapply(dir_names, function(dir) paste(dir, "/logs", sep=""))

HPD_list=list()
ESS_list=list()
n_good_runs=list()

for (dir in 1:length(dir_names)){
  
  
  
  
  repeat{
    
    dir_name = log_dirs[[dir]]
    
    list.files(paste(dir_name,"/" ,sep=""), all.files=TRUE)
    file_list = list.files(paste(dir_name,"/" ,sep=""))
    
    #prior_file_list=list.files(prior_dir)
    
    file_list=file_list[grep("_log$",unlist(file_list),fixed=FALSE)]
    
    #if using RFBS_DEC comp script
    run_list=file_list[grep("RFBS",unlist(file_list),fixed=FALSE)]
    
    run_list=run_list[-1]
    
    run_list = run_list[as.numeric(unlist(lapply(str_split(run_list, "_"), function(i) i[[1]])))<90]
    
    #run_list=file_list[grep("DEC",unlist(file_list),fixed=FALSE)]
    
    # prior_run=prior_file_list[grep("_log$",unlist(prior_file_list),fixed=FALSE)]
    #prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
    
    
    #sim_pars_strings = gsub("__","_",unlist(run_list),fixed=FALSE)
    
    sim_pars_strings = gsub("_log","",gsub("__","_",unlist(run_list),fixed=FALSE),fixed=FALSE)
    
    sim_pars_strings_unique=unique(sim_pars_strings)
    
    
    #sim_pars=as.numeric(unlist(str_split(sim_pars_strings_unique, "_")))
    
    sim_pars_full=do.call(rbind,lapply(sim_pars_strings_unique, function(run)  as.numeric(unlist(str_split(run, "_")))))
    
    if(is.null(sim_pars_full)==F){
      
      break
    }
  }
  
  sim_pars=sim_pars_full[,4:ncol(sim_pars_full)]
  
  sim_pars_off = sim_pars==0
  
  rowSums(sim_pars_off)
  
  burnin=0.5
  sim_pars_filler= sim_pars
  sim_pars_filler = sim_pars_filler * Inf
  {
    isin_HPD   =sim_pars_filler 
    lower_HPD  =sim_pars_filler 
    upper_HPD  =sim_pars_filler 
    ESS        =sim_pars_filler 
    post_mean  =sim_pars_filler 
    post_median =sim_pars_filler 
    switch_on_PD = sim_pars_filler 
    correct_switch_ratio = sim_pars_filler 
    switch_on = sim_pars !=0
    
    savage_dickey_ratios = sim_pars_filler 
    
    for (sim in 1:length(  sim_pars_strings_unique)){
      
      #  for (sim in 1:5){
      print(sim) 
      #sort files based on same simulating pars
      runs=run_list[sim_pars_strings==sim_pars_strings_unique[[sim]] ]
      
      {
      #col.names = c("iter",
      #              "lL_c",
      #              "prunelL_c",
      #              "rate_priorlL_c",
      #              "clado_priorlL_c",
      #              "_rf1_l_s",
      #             "_rf2_l_d",
      #              "_rf1_g_s",
      #              "_rf2_l_s",
      #              "_rf2_g_d",
      #              "_rf2_g_s",
      #              "switch_rf1_l_s",
      #             "switch_rf2_l_d",
      #              "switch_rf1_g_s",
      #              "switch_rf2_l_s",
      #              "switch_rf2_g_d",
      #              "switch_rf2_g_s",
      #              "split",
      #              "sub",
      #              "equal",
      #              "switch_split",
      #              "switch_sub",
      #              "switch_equal",
      #              "AR_rf1_l_s",
      #             "AR_rf2_l_d",
      #              "AR_rf1_g_s",
      #              "AR_rf2_l_s",
      #              "AR_rf2_g_d",
      #              "AR_rf2_g_s",
      #              "AR_switch_rf1_l_s",
      #             "AR_switch_rf2_l_d",
      #              "AR_switch_rf1_g_s",
      #              "AR_switch_rf2_l_s",
      #              "AR_switch_rf2_g_d",
      #              "AR_switch_rf2_g_s",
      #              "AR_clado_split_prob",
      #              "AR_clado_sub_prob",
      #              "AR_clado_equal_prob",
      #             "AR_switch_clado_split_prob",
      #             "AR_switch_clado_sub_prob",
      #             "AR_switch_clado_equal_prob")
      }
      
      
      #col.names = col.names[!grepl("l_d",col.names)]
      
      test=lapply(runs, function(run) read.table(paste(dir_name,run, sep="/"), header = T))
      
      for (run in 1:length(runs)){
        
        chain_length=nrow(test[[run]])
        
        test[[run]] = test[[run]][(chain_length*burnin):chain_length,]
        
      }
      
      

      if(chain_length<100){
        
        lower_HPD[sim,]=NA
        upper_HPD[sim,]=NA
        
        ESS[sim,]=NA
        isin_HPD[sim,]=NA
        post_mean[sim,]=NA
        post_median[sim,]=NA
        
        switch_on_PD[sim,] = NA
        correct_switch_ratio[sim,] = NA
        #ESS_list[[dir]]=ESS
        #HPD_list[[dir]]=isin_HPD
        
        
        
      }else{
        
        #plot(density(test[[1]]$sub))
        #polygon(density(test[[1]]$split))
        ##polygon(density(test[[1]]$sub))
        #polygon(density(test[[1]]$equal))
        
        npars= length(grep("X_" ,colnames(test[[run]]))) +3
        nswitches = npars
        #-4 without clado par,-5 with
        
        npars=(ncol(test[[1]])-5)/2
        
        #5:(4+) without clado par, 6:(5+)
        
        par_ind_vec   = c(grep("X_" ,colnames(test[[run]])), grep("^split$" ,colnames(test[[run]])),grep("^sub$" ,colnames(test[[run]])), grep("^equal$" ,colnames(test[[run]])))
        switch_ind_vec =grep("^switch_" ,colnames(test[[run]])) 
        
        chain_list=colnames(test[[1]])[par_ind_vec]
        
        
        
        
        for (chain in 1:length(chain_list)){
          #print(chain)
          HPD=HPDinterval(as.mcmc(test[[run]][,par_ind_vec[[chain]]]), prob=0.95)   
          
          #kd = density(test[[run]][,par_ind_vec[[chain]]] )
          #plot(kd, main=chain, ylab="yy")
          #
          #estimated_density <- approx(kd$x, kd$y, xout = 0)$y
          #estimated_density[is.na(estimated_density)] = 0.0
          
          y= test[[run]][,par_ind_vec[[chain]]] 
          xx <- seq(0, max(y), length = 500)
          ## bandwidth
          h <- 0.01
          ## get KDE using Gamma kernel
          den <- Gamma(x = xx, y = y, k = 500, h = h)
         # plot(den)
          ## evaluate at x=0
          estimated_density = den$y[den$x==0]
        
          
          if(par_ind_vec[[chain]]>17){
            prior_density = dunif(x = 0, min = 0, max = 1)
            
          }else{
            prior_density = dexp(x = 0, rate = 1)
          }
          
          savage_dickey_ratios[sim,chain] = log(estimated_density/prior_density)
          
          lower_HPD[sim,chain]=HPD[1]
          upper_HPD[sim,chain]=HPD[2]
          
          ESS[sim,chain]=effectiveSize(as.mcmc(test[[run]][,par_ind_vec[[chain]]]))   
          isin_HPD[sim,chain]=between(sim_pars[sim,chain],HPD[1],HPD[2])
          post_mean[sim,chain]=mean(test[[run]][,par_ind_vec[[chain]]])
          post_median[sim,chain]=median(test[[run]][,par_ind_vec[[chain]]])
          
          ESS_list[[dir]]=ESS
          HPD_list[[dir]]=isin_HPD
          
          
        }
        
        for (chain in 1:length(switch_ind_vec)){
          #print(chain)
          
          #sum(test[[run]][[,switch_ind_vec[[chain]]]]])
          switch_on[sim,chain] = sim_pars[sim,chain] != 0
          
          switch_on_PD[sim,chain] = sum(test[[run]][,par_ind_vec[[chain]]]!=0)/nrow(test[[run]])

          #test[[run]]$X_rf1_l_s!=0
          
          #test[[run]]$switch_rf1_l_s
          
        
          if(switch_on_PD[sim,chain]>1.0){
            print(sim)
            print(chain)
            XXXXXXX
          }
          
         if(sim_pars[sim,chain] == 0){
           
           correct_switch_ratio[sim,chain]  =   switch_on_PD[sim,chain]
           
         }else {
           
           correct_switch_ratio[sim,chain]  =  1 - switch_on_PD[sim,chain]
           
         }
          
          
        }
        
        
      }
      
      
      if(any(na.omit(correct_switch_ratio[sim,]) >1.0)){
        print(sim)
        print("632")
        print(chain)
        XXXXXXX
      }
      
      
    }
    
    
  }  
  
  if(length(which(correct_switch_ratio>1))>0){
    print(dir)
    print(which(correct_switch_ratio>1))
    print("641")
    XXXXXX
  }
  
  
  
  #mean(correct_switch_ratio[,1][sim_pars[,1]==0])
  
  #print(round(rbind(isin_HPD,colSums(isin_HPD)/nrow(isin_HPD)),digits = 4))
  
  #print(colSums(isin_HPD)/nrow(isin_HPD))
  print(dir_name)
  
  
  
  ESS_total=ESS
  HPD_total=HPD
  isin_HPD_total=isin_HPD
  
  
  
  sum(ESS_total<100)/length(ESS_total)
  
  post_mean_total=post_mean
  post_median_total=post_median
  sim_pars_total=sim_pars
  
  HPD_colors=c("red", "blue")
  
  

  sim_pars_total_list[[dir]]       =    sim_pars_total   
  post_median_total_list[[dir]]    =   post_median_total
  lower_HPD_list[[dir]]            =   lower_HPD       
  upper_HPD_list[[dir]]            =   upper_HPD       
  isin_HPD_list[[dir]]             =   isin_HPD        
  chain_list_list[[dir]]           =   chain_list      
  switch_on_PD_list[[dir]]         = switch_on_PD
  correct_switch_ratio_list[[dir]] = correct_switch_ratio
  switch_on_list[[dir]]            = switch_on
  savage_dickey_ratios_list[[dir]] = savage_dickey_ratios
  
  if(length(which(correct_switch_ratio_list[[dir]]>1))>0){
    print(dir)
    print(which(correct_switch_ratio_list[[dir]]>1))
    XXXXXX
  }
  
  
}


#par_names    = col.names[par_ind_vec ][1:5]
#switch_names = col.names[switch_ind_vec][1:5]

#####plot#####


correct_switch_ratio_list[[1]][12]

correct_switch_ratio_list[[1]][3,]
which(correct_switch_ratio_list[[1]]>0.8 & !switch_on_list[[1]])
sum(correct_switch_ratio_list[[1]]>1)
  
library(vioplot)

{
    

    #palette <- distinctColorPalette(n)
    
    palette=c("red", "blue", "gray" )
    

    col_vec=c("red", "blue")
    

    par_names=c(
     "1-0",
     "2-0", 
     "0-1",
     "2-2",
     "2-1",
     "0-2",
     "1-2"
    )
    
    for( i in 1:length(par_names)){
      
      pdf(paste("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim_plots/",par_names[[i]] ,"_ana_switch_support_test.pdf", sep=""),width = 8,height = 10)
      
      par(mfrow=c(3,2))
      
      for(x in 1:length(correct_switch_ratio_list)){


        
          
           on_correct    = (switch_on_PD_list[[x]][switch_on_list[[x]][,i]==1,i])
          off_correct    = (switch_on_PD_list[[x]][switch_on_list[[x]][,i]!=1,i])
         # either_correct = (correct_switch_ratio_list[[x]][,i])
       
          vioplot(on_correct,off_correct,names = c("on", "off"),col = col_vec,ylim = c(0,1), main= new_expr <- bquote(.(par_names[[i]]) ~ .(dir_labels[[x]])))

       }
        
      dev.off()
    }
}


library(vioplot)

{
  
  
  #palette <- distinctColorPalette(n)
  
  palette=c("red", "blue", "gray" )
  
  
  col_vec=c("red", "blue")
  
  
  par_names=c(
    "1-0",
    #"2-0", 
    "0-1",
    "2-2",
    "2-1",
    "0-2",
    "1-2"
  )
  
  for( i in 1:length(par_names)){
    
    pdf(paste("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim_plots/","sw_nodouble_",par_names[[i]] ,"_savage_dickey.pdf", sep=""),width = 8,height = 10)
    
    par(mfrow=c(3,2))
    
    max_lim=0
    min_lim=0
    for(x in 1:length(correct_switch_ratio_list)){
      on_correct    = (savage_dickey_ratios_list[[x]][switch_on_list[[x]][,i]==1,i])
      off_correct    = (savage_dickey_ratios_list[[x]][switch_on_list[[x]][,i]!=1,i])
      on_correct  = na.omit(on_correct[on_correct <Inf])
      off_correct = na.omit(off_correct[off_correct <Inf])
      
      
     max_lim =  max(c(max_lim, on_correct , off_correct ))
     min_lim =  min(c(max_lim, on_correct , off_correct ))
     
     max_lim =  3
     min_lim =  -3
     
     
    }
    for(x in 1:length(correct_switch_ratio_list)){
        
      on_correct    = (savage_dickey_ratios_list[[x]][switch_on_list[[x]][,i]==1,i])
      off_correct    = (savage_dickey_ratios_list[[x]][switch_on_list[[x]][,i]!=1,i])
      # either_correct = (correct_switch_ratio_list[[x]][,i])
      on_correct  = na.omit(on_correct[on_correct <Inf])
      off_correct = na.omit(off_correct[off_correct <Inf])
      vioplot(on_correct,off_correct,names = c(paste(length(on_correct), "on"), paste(length(off_correct), "off")),col = col_vec,ylim = c(min_lim,max_lim), main= new_expr <- bquote(.(par_names[[i]]) ~ .(dir_labels[[x]])))
      
    }
    
    dev.off()
  }
}



#post_median_total/




for( i in 1:length(par_names)){
  
  
  

  on_correct    = (switch_on_PD_list[[x]][switch_on_list[[x]][,i]==1,i])
  off_correct    = (switch_on_PD_list[[x]][switch_on_list[[x]][,i]!=1,i])
  # either_correct = (correct_switch_ratio_list[[x]][,i])
    
}



{
  par_names=c(
    expression(italic(l)[1 %->% 0]), 
    #expression(italic(l)[2 %->% 0]),
    expression(italic(g)[0 %->% 1]),
    expression(italic(sw)[2]),
    expression(italic(l)[2 %->% 1]),
    expression(italic(g)[0 %->% 2]),
    expression(italic(g)[1 %->% 2]),
    expression(italic(b)[i]),
    expression(italic(w)[i]),
    expression(italic(e)[i])
  )
  n_transpars = length( par_names) - 3
  
  ESS_rows=1:n_transpars 
  
  
  
  
  
 # pdf(paste(workdir,"/", "95p_3_150tip.pdf", sep=""), width = 20, height = 15)
  #par(mfrow=c(5,2))
  #layout_matrix= createLayoutMatrix(Nrow = 2, Ncol = 4, blockRows = 2, blockCols = 2, gridFillOrder  = "byrow", blockFillOrder = "byrow")     
  #layout(layout_matrix)
  
  # par(oma=c(0,0,3,0));  
  
  for (dir in 1:length(sim_pars_total_list)){
    
    pdf(paste(workdir,"/", dir, "_95p_3_150tip_nodouble_.pdf", sep=""), width = 20, height = 15)
    par(mfrow=c(5,2))
    
    # pdf(paste(workdir,"/",dir_name,"/",dir_name, "_posterior_median_lineplot.pdf", sep=""))
    
    
    
    sim_pars_total   =sim_pars_total_list[[dir]]   
    post_median_total=post_median_total_list[[dir]]
    lower_HPD       =lower_HPD_list[[dir]]        
    upper_HPD       =upper_HPD_list[[dir]]        
    isin_HPD        =isin_HPD_list[[dir]]         
    chain_list      =chain_list_list[[dir]]       
    ESS             = ESS_list[[dir]]
    
    
    
    
    
    #convert NA rows to 0 to filter them out with the ESS check, rather than writting a seperate filtering step
    ESS[is.na(ESS)]=0
    
    goodESS_rows=(1:nrow(ESS))[rowSums(do.call(cbind, (lapply(ESS_rows, function(t) ESS[,t]>200)))) == max(rowSums(do.call(cbind, (lapply(ESS_rows, function(t) ESS[,t]>200)))))]
    
    if(length(goodESS_rows)>=100){
      
      good_rows=sample(goodESS_rows,size = 100, replace = F)
      
    }else{
      
      good_rows=goodESS_rows
    }
    
    n_good_runs[[dir]]=length(good_rows)
    
    {
      #plot
      ncol_plots=4
      nrow_plots=2
      #par(mfrow=c( nrow_plots, ncol_plots ))
      
      
      #par(pty = "s")
      
      for (i in 1:ncol(post_median_total)){
        
        if(i>n_transpars){
          
          max_val= 1
          
          max=max_val+max_val*0.1
          maxy=max
          
          plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n", xlab = "", ylab="", xaxp = round(c(0, max_val, 2),digits = 0), yaxp = round(c(0, max_val, 2),digits = 0),  cex.main=2.5, cex.axis=2.5)
          
          
        }else{
          
          max_val=max(unlist(na.omit(post_median_total[post_median_total<Inf])))
          
          max=as.integer(max_val+max_val*0.1)
          
          if(max%%2==1){
            
            max=max+1
          }
          
          maxy=max
          
          plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n", xlab = "", ylab="", xaxp = as.integer((c(0, max, 2))), yaxp = as.integer((c(0, max, 2))),  cex.main=2.5, cex.axis=2.5)
          
          
        }  
        
        
        HPD_perc_round=round((colSums(isin_HPD[good_rows,])/nrow(isin_HPD[good_rows,]))[i],digits = 2 )
        HPD_perc=(colSums(isin_HPD[good_rows,])/nrow(isin_HPD[good_rows,]))[i]
        
        HPD_good_range=calc_HPD_ci(length( good_rows), 0.95)$ci
        if(HPD_perc>HPD_good_range[[1]] &  HPD_perc<HPD_good_range[[2]]){
          
          text(0+(max*.22), maxy-(maxy*.05), round((colSums(isin_HPD[good_rows,])/nrow(isin_HPD[good_rows,]))[i],digits = 2 ), col="blue", cex=3.5)
          
          # text(0+(max*.3), maxy-(maxy*.1), paste("n=", length(good_rows)), col="green", cex=3)
          
        } else{
          
          text(0+(max*.22), maxy-(maxy*.05), round((colSums(isin_HPD[good_rows,])/nrow(isin_HPD[good_rows,]))[i],digits = 2 ), col=2, cex=3.5)
          #text(0+(max*.3), maxy-(maxy*.1), paste("n=", length(good_rows)), col="green", cex=3)
          
          
        }
        points(sim_pars_total[good_rows,i], post_median_total[good_rows,i])
        
        
        
        arrows(sim_pars_total[good_rows,i], lower_HPD[good_rows,i], sim_pars_total[good_rows,i], upper_HPD[good_rows,i], length=0.01, angle=90, code=3,col =HPD_colors[isin_HPD[good_rows,i]+1] )
        #(ESS[,6]>200)
        #plot(sim_pars_total[,i], post_median_total[,i])
        abline(a=0,b=1)
      }
      
      
    }
    # post_median_plots[[dir]] <- recordPlot()
    
    # mtext(paste(dir_labels[[dir]]), side=3, line=0, cex=2 , outer=TRUE)  #"n=", length(good_rows))
    #box("outer", cex=2) 
    dev.off()
  }
  #file.copy(paste(workdir,"/",dir_name,"/", dir_name, "_posterior_median_lineplot.pdf", sep=""), "/Volumes/michael.landis/Active/RFBS_RIS/rf_sim_plots/posterior_medians/")
  #title(  paste(uncertainty_labels[[dir]]), line = -27, outer = TRUE)
  
  
  
  
  
 
  
}


#######

ESS_total=rbind(ESS_total,ESS)

HPD_total=rbind(HPD_total,HPD)
isin_HPD_total=rbind(isin_HPD_total,isin_HPD)
post_mean_total=rbind(post_mean_total,post_mean)
post_median_total=rbind(post_median_total,post_median)
sim_pars_total=rbind(sim_pars_total,sim_pars)
colSums(isin_HPD_total)/nrow(isin_HPD_total)


plot(pbtree(n = 3), )






#if((i-1)%%ncol_plots==0){
#  
#  if((i)>(nrow_plots*ncol_plots-ncol_plots)){
#    
#    plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n",xlab = "True Simulating Pars", ylab="Posterior Median", cex.main=2.5)
#    
#  }else{
#    
#    plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n", xlab="", ylab="Posterior Median",cex.main=2.5)
#    
#  }
#  
#}else{
#  
#  if(i>(nrow_plots*ncol_plots-ncol_plots)){
#    
#    plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n",xlab = "True Simulating Pars", ylab="", cex.main=2.5)
#    
#  }else{
#    
#    plot(c(0, max), c(0, maxy), main=par_names[[i]],type = "n", xlab = "", ylab="", cex.main=2.5)
#    
#  }
#  


#}

####messed up file naming when adding clado par, use as below
#
#file_list=list.files(dir_name)
#prior_file_list=list.files(prior_dir)
#
#run_list=file_list[grep("_log.txt",unlist(file_list),fixed=FALSE)]
#
#prior_run=prior_file_list[grep("__log.txt",unlist(prior_file_list),fixed=FALSE)]
#prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
#
#
#sim_pars_string_ = gsub("_log.txt*","",unlist(run_list),fixed=FALSE)
#sim_pars_strings = gsub("\t","",unlist(sim_pars_string_),fixed=FALSE)
#
#sim_pars_strings_unique=unique(sim_pars_strings)
#
#
#sim_pars=as.numeric(unlist(str_split(sim_pars_strings_unique, "_")))
#
#sim_pars=do.call(rbind,lapply(sim_pars_strings_unique, function(run)  as.numeric(unlist(str_split(run, "_")))))
#

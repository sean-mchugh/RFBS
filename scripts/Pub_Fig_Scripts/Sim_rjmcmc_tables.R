
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
  
  dir_names= c("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_rf_gl_ds")
               #'/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_unce[1.0][0.33][0.25]_rf_gl_ds',
               #'/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_unce[1.0][0.33][0.75]_rf_gl_ds',
               #'/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_unce[1.0][0.66][0.25]_rf_gl_ds',
               #'/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_1_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_unce[1.0][0.66][0.75]_rf_gl_ds')

  dir_names= c("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/sim/resub/rj_switch/BSim_RFBS_DEC_comp150t_3nB_Exp1p0iter_1000000_Ncladooff_0_Nrateoff_2_f2n_clado_Rtre_2g_2l_1g_1l_2sw_dr_rf_gl_ds")

    
  dir_labels=c("Full data",
               #"50 Tips / 1.0 rate prior",
               "0.33 biomes / 0.25 tips",
               "0.33 biomes / 0.75 tips",
               #"150 Tips / 1.0 rate prior",
               #"150 Tips / 5.0 rate prior",
               "0.66 biomes / 0.25 tips",
               "0.66 biomes / 0.75 tips"
  )
  
  
  #dir_labels=c("Full data")
  
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
  
  
  #plot_prior=T
  
  #prior_dir="BSim_RFBS_DEC_compPO_50t_3nB_Exp0p5iter_500000_f2n_clado_Rtre_2g_1g_1l_dr_rf_gl_ds"
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
mispsecified_on_list  = list()
mispsecified_off_list = list()

#savage_dickey_ratios_list=list()


log_dirs     =  lapply(dir_names, function(dir) paste(dir, "/logs", sep=""))

HPD_list=list()
ESS_list=list()
n_good_runs=list()



for (dir in 1:length(dir_names[[1]])){
  
  
  
  
  repeat{
    
    dir_name = log_dirs[[dir]]
    
    list.files(paste(dir_name,"/" ,sep=""), all.files=TRUE)
    file_list = list.files(paste(dir_name,"/" ,sep=""))
    
    #prior_file_list=list.files(prior_dir)
    
    file_list=file_list[grep("_log$",unlist(file_list),fixed=FALSE)]
    
    #if using RFBS_DEC comp script
    run_list=file_list[grep("RFBS",unlist(file_list),fixed=FALSE)]
    
    #run_list=run_list[-1]
    
    run_list = run_list[as.numeric(unlist(lapply(str_split(run_list, "_"), function(i) i[[1]])))]
    
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
    mispsecified_on = sim_pars_filler 
    mispsecified_off = sim_pars_filler 
    
    switch_on = sim_pars !=0
    
    savage_dickey_ratios = sim_pars_filler 
    
    for (sim in 1:length(  sim_pars_strings_unique)){
      
      #  for (sim in 1:5){
      print(sim) 
      #sort files based on same simulating pars
      runs=run_list[sim_pars_strings==sim_pars_strings_unique[[sim]] ]
      
      #col.names = col.names[!grepl("l_d",col.names)]
      
      test=lapply(runs, function(run) read.table(paste(dir_name,run, sep="/"), header = T))
      
      for (run in 1:length(runs)){
        
        chain_length=nrow(test[[run]])
        
        #test[[run]] = test[[run]][(chain_length*burnin):chain_length,]
        test[[run]] = test[[run]][500:1000,]
        
      }
      
      
      
      if(chain_length<1000){
        
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
          #den <- Gamma(x = xx, y = y, k = 500, h = h)
          # plot(den)
          ## evaluate at x=0
          #estimated_density = den$y[den$x==0]
          
          
          #if(par_ind_vec[[chain]]>17){
          #  prior_density = dunif(x = 0, min = 0, max = 1)
          #  
          #}else{
          #  prior_density = dexp(x = 0, rate = 1)
          #}
          
         # savage_dickey_ratios[sim,chain] = log(estimated_density/prior_density)
          
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
            mispsecified_off[sim,chain]  = NA
            
            if(correct_switch_ratio[sim,chain] < 0.05){
              mispsecified_on[sim,chain]  = TRUE
            }
          }else {
            
            correct_switch_ratio[sim,chain]  =  1 - switch_on_PD[sim,chain]
            mispsecified_on[sim,chain]  = NA
            if(correct_switch_ratio[sim,chain] < 0.05){
              mispsecified_off[sim,chain]  = TRUE
            }
            
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
  mispsecified_on_list[[dir]]  =  mispsecified_on
  mispsecified_off_list[[dir]] = mispsecified_off
  
 # savage_dickey_ratios_list[[dir]] = savage_dickey_ratios
  
  if(length(which(correct_switch_ratio_list[[dir]]>1))>0){
    print(dir)
    print(which(correct_switch_ratio_list[[dir]]>1))
    XXXXXX
  }
  
  
}


correct_switch_ratio_list[[1]][]


par_names=c(
  expression(paste("Enabled loss event "      , italic(l)[1 %->% 0])), 
  expression(paste("Dobule loss event "       , italic(l)[2 %->% 0])), 
  expression(paste("Enabled gain event "      , italic(g)[0 %->% 1])),
  expression(paste("Established switch event ", italic(sw)[2])),
  expression(paste("Established loss event "  , italic(l)[2 %->% 1])),
  expression(paste("Double gain event "       , italic(g)[0 %->% 2])),
  expression(paste("Established gain event "  , italic(g)[1 %->% 2]))
)



rj_on_PD = switch_on_PD_list[[1]][,6]
on       = switch_on_list[[1]][,6]
sig_thresh = 0.95

rj_sim_sum = function(rj_on_PD, on, sig_tresh){ 
  sig_on       = rj_on_PD[on]  > sig_thresh
  miss_sig_off = rj_on_PD[on]  < (1-sig_thresh)
  sig_off      = rj_on_PD[!on] < (1-sig_thresh)
  miss_sig_on  = rj_on_PD[!on] > (sig_thresh)
  
  # True ON
  true_on = sum(na.omit(sig_on ))/length(na.omit(sig_on ))
  
  # False On
  false_off = sum(na.omit(miss_sig_off ))/length(na.omit(miss_sig_off ))
  
  # True Off
  true_off = sum(na.omit(sig_off ))/length(na.omit(sig_off ))
  
  # False ON
  false_on = sum(na.omit(miss_sig_on ))/length(na.omit(miss_sig_on ))
  
  return(list(true_on=true_on, false_off=false_off ,true_off=true_off, false_on = false_on ))
}


rj_sim_sum <- function(rj_on_PD,
                       on,                # logical: TRUE  = actual positive, FALSE = actual negative
                       sig_thresh = 0.50  # probability threshold that decides ON vs OFF
){
  ## ------------------------------------------------------------------
  ## 1. Generate binary predictions
  ## ------------------------------------------------------------------
  pred_on <- rj_on_PD >= sig_thresh  # predicted positives
  pred_off <- !pred_on               # predicted negatives
  
  ## ------------------------------------------------------------------
  ## 2. Confusion-matrix cell counts
  ## ------------------------------------------------------------------
  TP <- sum(pred_on  &  on , na.rm = TRUE)
  FP <- sum(pred_on  & !on , na.rm = TRUE)
  FN <- sum(pred_off &  on , na.rm = TRUE)
  TN <- sum(pred_off & !on , na.rm = TRUE)
  
  ## Totals -----------------------------------------------------------
  P <- TP + FN                         # actual positives
  N <- TN + FP                         # actual negatives
  Tot <- P + N
  
  ## Helper to avoid division-by-zero warnings ------------------------
  div <- function(a, b) ifelse(b == 0, NA_real_, a / b)
  
  ## ------------------------------------------------------------------
  ## 3. Derived statistics (matching the chart)
  ## ------------------------------------------------------------------
  
  ## Basic rates
  prevalence <- div(P , Tot)                     # P  / (P+N)
  TPR <- div(TP , P )        # aka recall, sensitivity
  FNR <- div(FN , P )        # 1 – TPR
  TNR <- div(TN , N )        # specificity
  FPR <- div(FP , N )        # 1 – TNR
  
  ## Predictive values
  PPV <- div(TP , TP + FP)   # precision
  NPV <- div(TN , TN + FN)
  FDR <- 1 - PPV
  FOR <- 1 - NPV
  
  ## Accuracy family
  ACC <- div(TP + TN , Tot)
  BA  <- div(TPR + TNR , 2)
  F1  <- div(2 * TP , 2 * TP + FP + FN)
  
  ## Correlation / association indices
  FM  <- sqrt(PPV * TPR)                               # Fowlkes–Mallows
  MCC <- div( (TP*TN - FP*FN),
              sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) ) )
  MK  <- PPV + NPV - 1                                 # Markedness (Δp)
  
  ## Likelihood ratios & bookmakers
  LR_pos <- div(TPR , FPR)                             # LR+
  LR_neg <- div(FNR , TNR)                             # LR–
  BM     <- TPR + TNR - 1                              # Bookmaker informedness
  
  ## Diagnostic odds & prevalence threshold
  DOR <- div(LR_pos , LR_neg)
  PT  <- ifelse((TPR - FPR) == 0, NA_real_,
                (sqrt(TPR * FPR) - FPR) / (TPR - FPR))
  
  ## Threat / Jaccard
  TS <- div(TP , TP + FN + FP)                         # threat score, CSI, Jaccard
  
  ## ------------------------------------------------------------------
  ## 4. Return everything in a tidy list
  ## ------------------------------------------------------------------
  list(
    ## Confusion-matrix counts
    TP = TP, FP = FP, TN = TN, FN = FN,
    
    ## Primary rates
    prevalence = prevalence,
    TPR = TPR, FNR = FNR,
    TNR = TNR, FPR = FPR,
    
    ## Predictive values
    PPV = PPV, NPV = NPV,
    FDR = FDR, FOR = FOR,
    
    ## Accuracy & balanced metrics
    ACC = ACC, BA = BA, F1 = F1,
    
    ## Association / correlation
    FM = FM, MCC = MCC, MK = MK,
    
    ## Likelihood / diagnostic
    LR_pos = LR_pos, LR_neg = LR_neg,
    BM = BM, PT = PT, DOR = DOR,
    
    ## Threat / Jaccard
    TS = TS
  )
}


sig_thresh  = 0.8
rj_sum_list = list()
rj_on_PD    = switch_on_PD_list[[1]][,1:7]
on          = switch_on_list[[1]][,1:7]

rj_sum_list[[1]] = rj_sim_sum (rj_on_PD, on, sig_thresh)

for (i in 1:7){
  on       = switch_on_list[[1]][,i]
  rj_on_PD = switch_on_PD_list[[1]][,i]
  rj_sum_list[[i+1]] = rj_sim_sum (rj_on_PD, on, sig_thresh)
}


summary_mat= do.call(rbind, lapply(rj_sum_list, function(i) unlist(i)))
rownames(summary_mat) = c("All Events",par_names)
write.csv(file = "~/Downloads/miss_test.csv", x=summary_mat )


summary_mat

library(gridExtra)
library(grid)
# Plot as a table with expressions
grid.table(summary_mat, rows = rownames(summary_mat) , cols=colnames(summary_mat), parse=T )
#par_names    = col.names[par_ind_vec ][1:5]
#switch_names = col.names[switch_ind_vec][1:5]

library(grid)





library(ggplot2)
library(tidyr)
library(dplyr)

## names as *character* plot-math strings (parse needs strings)
par_names = c(  "Enabled loss",
  "Double loss",
  "Enabled gain",
  "Established switch",
  "Established loss",
  "Double gain",
  "Establishedgain"
)

col_math <- c(
  "italic(p)[true~on]",
  "italic(p)[true~off]",
  "italic(p)[false~on]",
  "italic(p)[false~off]"
)

x=summary_mat

## long form for ggplot
tbl <- as.data.frame(x) |>
  mutate(row = factor(seq_along(row_math), levels = rev(seq_along(row_math))),
         row_lbl = row_math) |>
  pivot_longer(cols = everything() & !row & !row_lbl,
               names_to = "col", values_to = "val") |>
  mutate(col = factor(col, levels = colnames(x)),
         col_lbl = col_math[as.numeric(col)])

ggplot(tbl) +
  geom_tile(aes(col, row), fill = NA, colour = "grey70") +
  geom_text(aes(col, row, label = sprintf("%.6f", val)), size = 3) +
  # row labels (x = 0 on an extra fake column)
  geom_text(aes(x = 0, y = row, label = row_lbl), parse = TRUE,
            hjust = 1, size = 3.5) +
  # column labels at top
  geom_text(data = distinct(tbl, col, col_lbl),
            aes(col, y = max(as.numeric(tbl$row)) + 1, label = col_lbl),
            parse = TRUE, fontface = "bold", size = 4) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void()




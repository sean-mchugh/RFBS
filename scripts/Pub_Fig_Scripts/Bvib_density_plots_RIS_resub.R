match_string_elements_to_string_vec <- function(string_vec, string_pieces, name_vec) {
  # Create an empty list to store the results
  scores_list <- vector("list", length(string_vec))
  
  # Loop over each directory name
  for (i in seq_along(string_vec)) {
    dir_name <- string_vec[i]
    # Check which evidence strings are present in the directory name
    scores_list[[i]] <- paste(name_vec[sapply(string_pieces, function(string_piece) grepl(string_piece, dir_name))], collapse=".")
  }
  
  return(scores_list)
}


library(stringr)
library(dplyr)
library(coda)
library(ape)
library(phytools)

library(randomcoloR) 
library(vioplot)

RFBS_dir_names=c("Bvib_3nB_LN_0_p5_1000000_admat_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_admat_excf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_admat_incf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_admat_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_admat_incf_excf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_excf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_incf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
                 "Bvib_3nB_LN_0_p5_1000000_incf_excf_eco_allo_clado_1g_1l_rf_gl",
                 "Bvib_3nB_LN_0_p5_1000000_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds")


RFBS_dir_names= 
  c("Bvib_3nB_LN_0_p5_1000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_admat_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
    "Bvib_3nB_LN_0_p5_1000000_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds")

RFBS_dir_names= 
  c( "Bvib_3nB_Exp0p5_10000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_admat_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
     "Bvib_3nB_Exp0p5_10000000_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds"
  )


RFBS_dir_names= 
  c(  "Bvib_3nB_Exp0p5_10000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_admat_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_incf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_10000000_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.conservative.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.conservative.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.conservative.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.germination.only.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.bold.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_excf_3.biomes.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds"
  )



#RFBS_dir_names= 
#  c(  "Bvib_3nB_Exp0p5_10000000_admat_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_10000000_admat_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.germination.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.germination.only.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.bold.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.conservative_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds",
#      "Bvib_3nB_Exp0p5_3000000_admat_excf_3.biomes.USDA_eco_allo_clado_2g_1g_1l_rf_gl_ds"
#  )

RFBS_dir_names = c(
  "Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa"
  
)

RFBS_dir_names = c(
  "Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_3000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa"
  
)


RFBS_dir_names = c(
  #"Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_admat_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_incf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_bold_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_excf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss",
  "Bvib_3nB_Exp0p5_1000000_Foss_noincf_cons_noexcf_3.biomes.germination.only.leafing.conservative.USDA_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa"
  
)


RFBS_dir_names = RFBS_dir_names[c(grep("RJa", RFBS_dir_names ))] 

#new_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_3000000",list.files()  ))]
#old_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_10000000",list.files()  ))][c( 5,1,6,7, 2,3,8,4)]

#RFBS_dir_names=c(old_RFBS_dir_names, new_RFBS_dir_names)

#"excf_3.biomes.germination.only.leafing.conservative.USDA"
#"excf_3.biomes.germination.only.germination.bold.leafing.bold"
#
#
#evidence_vec=c("germination.only" ,"germination.conservative" ,"germination.bold", "leafing.conservative", "leafing.bold","USDA", "admat", "_incf_bold", "_incf_cons","_excf"  )
#name_vec=c("G_O" ,"G_C" ,"G_B", "L_C", "L_B","USDA", "A", "Ib", "Ic", "E" )


evidence_vec=c( "Bvib_3nB_Exp0p5_1000000_PO_Foss_admat_incf_bold_excf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_2l_1g_1l_2sw_rf_gl_ds_Foss_RJa",
                "admat", "_incf_bold", "_incf_cons", "_excf_3.biomes.germination.only.germination.bold.leafing.bold", "_excf_3.biomes.germination.only.leafing.conservative.USDA"  )
name_vec=c("Prior", "A", "Ib", "Ic", "Eb", "Ec" )


RFBS_dir_names_sort=sort(RFBS_dir_names)



process_string <- function(x) {
  x <- gsub("admat", "A", x)  
  x <- gsub("excf", "E", x)   
  x <- gsub("incf", "I", x)   
  x <- gsub("_", "", x)   
  return(x)
}


# Apply the custom function to each element in the list using sapply or lapply
#RFBS_legend_names <- sapply(RFBS_legend_names, process_string)
#
#RFBS_legend_names[is.na(RFBS_legend_names)]="None"
#RFBS_legend_names[[length(RFBS_legend_names)+1]]="Prior"

#RFBS_legend_order=c(9, 5,1,6,7, 2,3,8,4)
#RFBS_legend_order=1:length(RFBS_legend_names)


#old_RFBS_legend_names <- lapply(RFBS_dir_names[1:9], function(dir) str_match(dir, "10000000_\\s*(.*?)\\s*_eco")[[2]])

#RFBS_legend_names <- lapply(RFBS_dir_names, function(dir) str_match(dir, "3000000_\\s*(.*?)\\s*_eco")[[2]])


#RFBS_legend_names <-lapply(RFBS_dir_names, function(dir) str_match(dir, "3000000_excf_3.biomes.\\s*(.*?)\\s*_eco")[[2]])


RFBS_legend_names = match_string_elements_to_string_vec(string_vec    = RFBS_dir_names_sort, 
                                                        string_pieces = evidence_vec  , 
                                                        name_vec)
#RFBS_legend_names[[33]] = "Prior" 

RFBS_legend_names[RFBS_legend_names == ""] = "None"

#RFBS_legend_order = c(9, 5,1,6,7, 2,3,8,4)

# Flatten the list into a character vector
RFBS_legend_names <- unlist(RFBS_legend_names)

RFBS_legend_names  = gsub(".", "_", RFBS_legend_names, fixed = T )

new_RFBS_runs_df = cbind(RFBS_legend_names, RFBS_dir_names_sort)

new_RFBS_runs_df = new_RFBS_runs_df[!duplicated(new_RFBS_runs_df[,1]), ]


#RFBS_treatments_sorted=c("Prior", "None", "A", "Ic","Ib", "Ec", "Eb", "A_Ic", "A_Ib", "A_Ec", "A_Eb","Ic_Ec", "Ic_Eb", "Ib_Ec", "Ib_Eb", "A_Ic_Ec", "A_Ic_Eb", "A_Ib_Ec", "A_Ib_Eb")
RFBS_treatments_sorted=c("None", "A", "Ic","Ib", "Ec", "Eb", "A_Ic", "A_Ib", "A_Ec", "A_Eb","Ic_Ec", "Ic_Eb", "Ib_Ec", "Ib_Eb", "A_Ic_Ec", "A_Ic_Eb", "A_Ib_Ec", "A_Ib_Eb")


sorted_RFBS_runs_df= new_RFBS_runs_df[match(RFBS_treatments_sorted,new_RFBS_runs_df[,1] ) , ]

RFBS_dir_names = sorted_RFBS_runs_df[-1,2]
RFBS_dir_names = sorted_RFBS_runs_df[ ,2]

RFBS_Prior_dir = "Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss"
  sorted_RFBS_runs_df[1,2]

RFBS_legend_names = new_RFBS_runs_df[,1]


#setwd("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp/viburnum/resub")

workdir="/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp/viburnum/resub/"

RFBS_dir_names = paste0(workdir, RFBS_dir_names)
#RFBS_Prior_dir = paste0(workdir, RFBS_Prior_dir)

#dir_name=dir_names[[1]]

tree=read.tree("/Volumes/michael.landis/Active/Sean/RFBS/data/emp/viburnum/viburnum_sorted.tre")

rescale=max(nodeHeights(tree))
tree$edge.length

######make rfbs posterior objects

sim_dat=F
post_dist=list()
i=0



  
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
  
  
  # dir_name=dir_names[[i]]
  
  RFBS_file_list=lapply(RFBS_dir_names, function(dir) list.files(dir))
  #DEC_file_list=lapply(DEC_dir_names, function(dir) list.files(dir))
  
  
  ###
  prior_file_list=list.files(RFBS_Prior_dir)
  
  prior_run=prior_file_list[grep("__log.txt",unlist(prior_file_list),fixed=FALSE)]
  ####
  
  #prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
  #
  RFBS_run_files=lapply(RFBS_file_list, function(file_list) file_list[grep("_log$",unlist(file_list),fixed=FALSE)])
  #RFBS_prior_run=prior_file_list[grep("_log$",unlist(prior_file_list),fixed=FALSE)]
  #RFBS_prior_test=read.table(paste(RFBS_Prior_dir,RFBS_prior_run[[1]], sep="/"), header = T)
  RFBS_prior_test = readRDS("prior_only_nodouble.RDS")
  burnin=0.5
  ####make rfbs post objects#####    
  ######this loads the chains takes awhile and need connection to RIS
  RFBS_chains_full=lapply(1:length(RFBS_dir_names), function(dir) lapply(RFBS_run_files[[dir]], function(run) read.table(paste(RFBS_dir_names[[dir]],run, sep="/"), header = T) ))
  ###############################################
  
  
{
  
  RFBS_chains=RFBS_chains_full
  #npars=19 #(ncol(RFBS_chains[[1]][[1]] )-4)/2
  
  npars=20 #(ncol(RFBS_chains[[1]][[1]] )-4)/2
  

  #chain_ind_vec=6:(5+npars)
  
  #chain_ind_vec = c(6:11, 18:20)
  chain_ind_vec = c(6:22)
  #chain_ind_vec = c(6:20)
  
  chain_list=colnames(RFBS_chains[[1]][[1]])[chain_ind_vec]
  RFBS_post_dist=list()
  
  {
    for (dir in (1:length(RFBS_chains))){
      print(dir)
      RFBS_post_dist[[dir]]=list()
      for (run in 1:length(RFBS_chains[[dir]])){
        chain_length=nrow(RFBS_chains_full[[dir]][[run]])
        RFBS_chains[[dir]][[run]] = RFBS_chains_full[[dir]][[run]][(chain_length*burnin):chain_length,]

      }
      for (chain in 1:length(chain_ind_vec)){
        RFBS_post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply( 1:length(RFBS_chains[[dir]]), function(run) RFBS_chains[[dir]][[run]][,chain_ind_vec[[chain]] ]))
      }
    }
  }
}

  
  
#chain_par_vec = c(6:11, 18:20)-5
#chain_par_names = chain_list[chain_par_vec]
#chain_rj_vec = c(12:17) - 5
  

chain_par_vec = c(6:12, 20:22)-5
chain_par_names = chain_list[chain_par_vec]
chain_rj_vec = c(13:19) - 5



RFBS_post_dist[[dir]][chain_rj_vec]

unique(as.numeric(RFBS_post_dist[[dir]][[chain]]))

RFBS_post_dens=list()
RFBS_rj_dens=list()
RFBS_names = names(RFBS_post_dist[[1]])

for( dir in (1:length(RFBS_post_dist))){
  RFBS_post_dens[[dir]]=list()
  RFBS_rj_dens[[dir]]=list()
  
  for( i in 1:length(chain_par_vec)){
    chain = chain_par_vec[[i]]
    RFBS_post_dens[[dir]][[i]]    =density(as.numeric(RFBS_post_dist[[dir]][[RFBS_names[chain]]]))
    names(RFBS_post_dens[[dir]][i])=   names(RFBS_post_dist[[dir]])[chain]
    
  }
  #names(RFBS_post_dens[[dir]])=   names(RFBS_post_dist[[dir]])[chain_par_vec]
  for(  i in 1:length(chain_rj_vec)){
    chain = chain_rj_vec[[i]]
    
    RFBS_rj_dens[[dir]][[i]]=density(as.numeric(RFBS_post_dist[[dir]][[RFBS_names[chain]]]))
    names(RFBS_rj_dens[[dir]][i])=   names(RFBS_post_dist[[dir]])[chain]
    
  }
}

prior_chains=list()

######plot RFBS########


RFBS_prior_test = readRDS("prior_only_nodouble.RDS")

{
  pdf(paste("~/Projects/RFBS-main/outfiles/emp/viburnum/resub_figs/test_RJa_conserv_doubleloss_I_3biome_bold_consv_EXP_0p5_bvib_RFBS_density_violin.pdf", sep=""),width = 10,height = 10)
  
  
  {
    {
      
      
      
      #par(mfrow=c(3,2))
      # par(oma = c(4,1,1,1), mfrow = c(4, 2), mar = c(2, 2, 2, 2))
      # Increase the bottom outer margin to provide more space for labels
      par(oma = c(6, 1, 1, 1))  # Increase the bottom outer margin
      
      # Increase the bottom margin of each plot to prevent label overlap
      # The 'mar' parameter takes the form c(bottom, left, top, right)
      par(mar = c(5, 2, 2, 2))  # Increase the bottom margin
      
      # Set the layout of the plotting area to 4x2
      par(mfrow = c(5, 2))
      

      par_names=c(
        expression(paste("Enabled loss event "      , italic(l)[1 %->% 0])), 
        expression(paste("Double loss event "       , italic(l)[2 %->% 0])), 
        expression(paste("Enabled gain event "      , italic(g)[0 %->% 1])),
        expression(paste("Established switch event ", italic(s)[2])),
        expression(paste("Established loss event "  , italic(l)[2 %->% 1])),
        expression(paste("Double gain event "       , italic(g)[0 %->% 2])),
        expression(paste("Established gain event "  , italic(g)[1 %->% 2])),

        expression(paste("Speciation event "      , italic(b))          ),
        expression(paste("Speciation event "      , italic(s))          ),
        expression(paste("Speciation event "      , italic(e))          )
      )
      
      RFBS_legend_names=c("Prior", "None", "A", "Ic","Ib", "Ec", "Eb", "A_Ic", "A_Ib", "A_Ec", "A_Eb","Ic_Ec", "Ic_Eb", "Ib_Ec", "Ib_Eb", "A_Ic_Ec", "A_Ic_Eb", "A_Ib_Ec", "A_Ib_Eb")
      
      # col_vec=c( "grey", "gray31", "blue4", "yellow3", "red1", "red4", "green3","purple1", "purple4","orange2", "darkorange2",  "burlywood3","chocolate4")
      
      
      col_vec=c( "ivory", "linen", "blue4", "yellow", "yellow3", "red1", "red4", "green1", "green4", "purple1", "purple4","orange1", "darkorange2", "orange3", "darkorange4",  "gray70", "gray50", "gray30", "gray0")
      
      prior_col="grey"
      
      #for (z in 1:9){
      for (z in 1:10){
          
        
        #old_par_ind = c(2,1,6,4,5,3,7,8, 9)[[z]]
        old_par_ind = c(3,1,7,5,6,2,4,8,9, 10)[[z]]
        
        chain=chain_par_vec[old_par_ind]
        
        #chain=c(10, 1:9)[[z]]
        
        

        #ymax= c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5,   1.0,1.0,1.0)
        ymax= c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.0,1.0,1.0)
        
        
        
        
        #plot(post_dist, main=chain_list[[chain]])
        if(old_par_ind <8){
          prior_chain = "X_rf1_l_s"
          
        }else{
          
          prior_chain = "equal"
        }
        prior_chain=unlist(sample(RFBS_prior_test[,prior_chain], nrow(RFBS_prior_test), replace = T))
        #prior_chain=readRDS("prior_no_dl.RDS")
        
        prior_dist=density(prior_chain)
        

        min_dist_size= min(unlist(c(lapply(1:length(RFBS_chains), function(dir) length(RFBS_post_dist[[dir]][[chain]])), length(unlist(prior_chain)))))
        #min_dist_size= min(unlist(c(lapply(1:length(RFBS_chains), function(dir) length(RFBS_post_dist[[dir]][[chain]])) )))
        
        #polygon(prior_dist, col=NA, border=t_col(prior_col,80) )
        dists=cbind(do.call(cbind, c(
          list(sample(unlist(prior_chain), size = min_dist_size, replace=F)),
          lapply(1:length(RFBS_chains), function(dir) sample(as.numeric(RFBS_post_dist[[dir]][[chain]]), size = min_dist_size, replace=F) ) ) )
        )
        
        which(is.na(dists), arr.ind=T)
        

        vioplot(dists,
                col=col_vec,
                names= rep("",19), #rep("",13),#rep("",9)
                main=par_names[[old_par_ind]],
                ylim=c(0,ymax[[z]]),
                
                #names=RFBS_legend_names,#rep("",30),#rep("",9)
                cex.main=2,
                cex.axis=1.0,
                #xlab=RFBS_legend_names[RFBS_legend_order] ,
                las=2
                
        )

      } 
      
    }
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")  

    legend('bottom',
           legend = RFBS_legend_names, 
           col = "black",  # Border color
           pt.bg = col_vec,  # Fill color for points
           pch = 21,  # Point character with border and fill
           pt.cex = 2,  # Point size
           cex = 2, 
           seg.len = 0.25, 
           ncol = 5, 
           bty = 'n')
  }
  
  dev.off()
}















{
  
  pdf(paste("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp_plots/viburnum/3correctIEA_rescaled_EXP_0p5_bvib_RFBS_DEC_density.pdf", sep=""),width = 8,height = 10)
  
  
  
  n <- 8
  palette <- distinctColorPalette(n)
  
  palette=c("green", "yellow", "cyan", "black", "grey", "red4", "blue4", "magenta" )
  
  xmax_vec=c(6.0,6.0,10.0,25,10.0,1.0)/rescale
  
  col_vec=c("yellow3", "orange", "green", "brown", "dark gray", "red4", "blue4", "purple" )
  
  prior_col="gray"
  
  col_vec[[length(col_vec)+1]]=prior_col
  
  
  {
    
    
    
    par(mfrow=c(3,2))
    
    
    par_names=c(
      "RFBS 1->0",
      "RFBS 0->1",
      "RFBS 2->1",
      "RFBS 0->2",
      "RFBS 1->2",
      "RFBS 2->2",
      "DEC and RFBS clado split prob"
      
    )
    
    for (z in 1:6){
      
      chain=c(2,1,5,3,4,6)[[z]]
      
      #plot(post_dist, main=chain_list[[chain]])
      prior_dist=density(unlist(sample(RFBS_prior_test[,chain_par_vec[[chain]]]/rescale, nrow(RFBS_prior_test), replace = T)))
      
      #abline(v=sim_pars[sim,chain],col=2)
      x_min=min(unlist(lapply(RFBS_post_dist, function(i) min( i[[chain]]))))
      x_min=x_min-x_min*.3
      x_max=max(unlist(lapply(RFBS_post_dist, function(i) max( i[[chain]]))))
      if(chain!=length(chain_par_vec)){
        x_max=x_max/5
        x_max=xmax_vec[[z]]
      }
      
      y_min=min(c(unlist(lapply(RFBS_post_dist, function(i) min( density(i[[chain]])$y   )))), min(  prior_dist$y)   )
      y_min=y_min-y_min*.3
      y_max=max(c(unlist(lapply(RFBS_post_dist, function(i) max( density(i[[chain]])$y   )))), max(  prior_dist$y)   )
      y_max=y_max+y_max*.3
      
      
      plot( c(0, x_max),c(0, y_max), main=par_names[[chain]],type = "n")
      
      
      
      
      #polygon(prior_dist, col=t_col(prior_col,50) )
      
      
      if(chain==2){
        
        #        legend("topright", legend=RFBS_legend_names, fill =unlist(lapply(col_vec, function(c) t_col(c, 20))))
        
      }
      
      
      for( dir in 1:length(RFBS_chains)){
        
        
        #polygon(density(RFBS_post_dist[[dir]][[chain]]), col=t_col(col_vec[[dir]],70))
        
        #if(chain==6){
        #  
        #  polygon(density(DEC_post_dist$clado_split_prob), border="orange", col=t_col("orange",20))
        #  
        #  legend("topright", legend=c(RFBS_legend_names,"DEC"), fill =c(unlist(lapply(col_vec, function(c) t_col(c, 20))), "orange"))
        #  
        #  
        #}
        
      }
      
      #polygon(prior_dist, col=NA, border=t_col(prior_col,80) )
      
      
      for( dir in 1:length(RFBS_chains)){
        
        
        polygon(density(RFBS_post_dist[[dir]][[chain]])
                , col=NA,lwd=0.5, border = t_col(col_vec[[dir]],0),  lwd=5 )
        
        
      }
      
      # polygon(prior_dist, col=NA, border=t_col(prior_col,0), lwd=5 )
      
      legend('topleft',legend=RFBS_legend_names, 
             col =unlist(lapply(col_vec, function(c) t_col(c, 20))),
             lwd = 5, xpd = TRUE, cex = 1, seg.len=1, bty = 'n', ncol=4)
      
      
    } 
    
    #par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    #plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    #plot_colors <- c("blue","black", "green", "orange", "pink")
    #legend(x = "top",inset = 0,
    #       legend = c("Fabricated Metal", "Iron and Steel", "Paper","Beverages", "Tobacco"), 
    #       col=plot_colors, lwd=5, cex=.5, horiz = TRUE)
    #legend(x = "top",inset = 0, legend=RFBS_legend_names, fill =unlist(lapply(col_vec, function(c) t_col(c, 20)))
    #        ,cex=.5, horiz = TRUE)
    # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    # plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    # legend('bottom',legend=RFBS_legend_names, 
    #        col =unlist(lapply(col_vec, function(c) t_col(c, 20))),
    #        lwd = 5, xpd = TRUE, cex = 1, seg.len=1, bty = 'n', ncol=4)
  }
  dev.off()
}





RFBS_post_dist=RFBS_post_dist 
RFBS_legend_names=RFBS_legend_names  
RFBS_legend_names[RFBS_legend_names == ""] = "None"

RFBS_chains = RFBS_chains


{
  
  
  
  pdf(paste("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp_plots/viburnum/3biome_all_treatments_Donoghue_EXP_0p5_bvib_RFBS_density_violin.pdf", sep=""),width = 10,height = 10)
  
  
  

  
  n <- 8
  palette <- distinctColorPalette(n)
  
  palette=c("green", "yellow", "cyan", "black", "grey", "red", "blue", "magenta" )
  
  xmax_vec=c(6.0,6.0,10.0,25,10.0,1.0)/rescale
  
  ymax=c(rep(1.5, 5),rep(1.0, 3))
  
  
  col_vec=c("yellow3", "orange2", "green4", "chocolate4", "gray31", "red4", "blue4", "purple4" )
  
  #col_vec=c("yellow", "orange", "green", "brown", "gray31", "red", "blue", "purple" )
  
  prior_col="gray"
  
  col_vec[[length(col_vec)+1]]=prior_col
  
  
  
  {
    {
      
      
      
      #par(mfrow=c(3,2))
      # par(oma = c(4,1,1,1), mfrow = c(4, 2), mar = c(2, 2, 2, 2))
      # Increase the bottom outer margin to provide more space for labels
      par(oma = c(6, 1, 1, 1))  # Increase the bottom outer margin
      
      # Increase the bottom margin of each plot to prevent label overlap
      # The 'mar' parameter takes the form c(bottom, left, top, right)
      par(mar = c(5, 2, 2, 2))  # Increase the bottom margin
      
      # Set the layout of the plotting area to 4x2
      par(mfrow = c(4, 2))
      
      par_names=c(
        expression(italic(l)[1 %->% 0]), 
        expression(italic(g)[0 %->% 1]),
        expression(italic(s)[2]),
        expression(italic(l)[2 %->% 1]),
        expression(italic(g)[0 %->% 2]),
        expression(italic(g)[1 %->% 2]),
        expression(italic(b)[i]),
        expression(italic(w)[i]),
        expression(italic(e)[i])
      )
      
      
      for (z in 1:9){
        
        chain=c(2,1,5,4,6,3, 7,8,9)[[z]]
        ymax= c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5,1.0,1.0,1.0)
        #plot(post_dist, main=chain_list[[chain]])
        prior_chain=unlist(sample(RFBS_prior_test[,chain_ind_vec[[chain]]], nrow(RFBS_prior_test), replace = T))
        prior_dist=density(prior_chain)
        
        #abline(v=sim_pars[sim,chain],col=2)
        #x_min=min(unlist(lapply(RFBS_post_dist, function(i) min( i[[chain]]))))
        #x_min=x_min-x_min*.3
        #x_max=max(unlist(lapply(RFBS_post_dist, function(i) max( i[[chain]]))))
        #if(chain!=length(chain_ind_vec)){
        #  x_max=x_max/5
        #  x_max=xmax_vec[[z]]
        #}
        
        y_min=min(c(unlist(lapply(RFBS_post_dist, function(i) min( density(as.numeric(i[[chain]]))$y   )))))#, min(  prior_dist$y)   )
        y_min=y_min-y_min*.3
        y_max=max(c(unlist(lapply(RFBS_post_dist, function(i) max( density(as.numeric(i[[chain]]))$y   )))))#, max(  prior_dist$y)   )
        y_max=y_max+y_max*.3
        
        #ymax= c(0.5, 0.5, 1.0, 1.5, 0.5, 1.0,1.0,1.0)
        
        
        #plot( c(0, x_max),c(0, y_max), main=par_names[[chain]],type = "n")
        
        
        
        
        #polygon(prior_dist, col=t_col(prior_col,50) )
        
        
        if(chain==2){
          
          #        legend("topright", legend=RFBS_legend_names, fill =unlist(lapply(col_vec, function(c) t_col(c, 20))))
          
        }
        
        
        for( dir in (1:length(RFBS_chains))){
          
          print(dir)
          
          #polygon(density(RFBS_post_dist[[dir]][[chain]]), col=t_col(col_vec[[dir]],70))
          
          #if(chain==6){
          #  
          #  polygon(density(DEC_post_dist$clado_split_prob), border="orange", col=t_col("orange",20))
          #  
          #  legend("topright", legend=c(RFBS_legend_names,"DEC"), fill =c(unlist(lapply(col_vec, function(c) t_col(c, 20))), "orange"))
          #  
          #  
          #}
          
        }
        
        #polygon(prior_dist, col=NA, border=t_col(prior_col,80) )
        dists=cbind(do.call(cbind, c(lapply(1:length(RFBS_chains), function(dir) as.numeric(RFBS_post_dist[[dir]][[chain]])) ) ),unlist(prior_chain))
        
        
        vioplot(dists,#[,RFBS_legend_order],
                col=col_vec, #[RFBS_legend_order],
                #names=RFBS_legend_names,
                main=par_names[[chain]],
                ylim=c(0,ymax[[z]]),
                names=RFBS_legend_names,#rep("",30),#rep("",9)
                cex.main=2,
                cex.axis=0.5,
                #xlab=RFBS_legend_names[RFBS_legend_order] ,
                las=2
                
        )
        
        #for( dir in 1:length(RFBS_chains)){
        #  
        #  
        #  polygon(density(RFBS_post_dist[[dir]][[chain]])
        #          , col=NA,lwd=0.5, border = t_col(col_vec[[dir]],0),  lwd=5 )
        #  
        #  
        #}
        
        #polygon(prior_dist, col=NA, border=t_col(prior_col,0), lwd=5 )
        
        #legend('topleft',legend=RFBS_legend_names, 
        #       col =unlist(lapply(col_vec, function(c) t_col(c, 20))),
        #       lwd = 5, xpd = TRUE, cex = 1, seg.len=1, bty = 'n', ncol=4)
        #
        
      } 
      
      #plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      #plot_colors <- c("blue","black", "green", "orange", "pink")
      #legend(x = "top",inset = 0,
      #       legend = c("Fabricated Metal", "Iron and Steel", "Paper","Beverages", "Tobacco"), 
      #       col=plot_colors, lwd=5, cex=.5, horiz = TRUE)
      #legend(x = "top",inset = 0, legend=RFBS_legend_names, fill =unlist(lapply(col_vec, function(c) t_col(c, 20)))
      #        ,cex=.5, horiz = TRUE)
      # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      # plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
      # legend('bottom',legend=RFBS_legend_names, 
      #        col =unlist(lapply(col_vec, function(c) t_col(c, 20))),
      #        lwd = 5, xpd = TRUE, cex = 1, seg.len=1, bty = 'n', ncol=4)
    }
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")  
    #legend('bottom',
    #     legend = RFBS_legend_names[RFBS_legend_order], 
    #     col = col_vec[RFBS_legend_order], lwd = 5, xpd = TRUE, cex = 2, seg.len=0.3,ncol = 9, bty = 'n')
  }
  
  dev.off()
}







######plot DEC########

DEC_chain_ind_vec


DEC_par_names=c("DEC 0->2", "DEC 2->0")

for (chain in c(1,2)){
  
  
  
  #plot(post_dist, main=chain_list[[chain]])
  #prior_dist=density(unlist(sample(DEC_prior_test[,chain_ind_vec[[chain]]], nrow(DEC_prior_test), replace = T)))
  
  #abline(v=sim_pars[sim,chain],col=2)
  x_min=min(DEC_post_dist[[chain]])
  x_min=x_min
  x_max=max((DEC_post_dist[[chain]]))
  x_max=x_max/1.6
  
  y_min=min(unlist(lapply(DEC_post_dist, function(i) min( density(i)$y   ))))  
  y_min=y_min-y_min*.3
  y_max=max(c(unlist(lapply(DEC_post_dist, function(i) max( density(i)$y   ))))  )
  y_max=y_max+y_max*.3
  
  
  plot( c(0, x_max),c(0, y_max), main=DEC_par_names[[ chain]],type = "n")
  
  
  
  
  
  #polygon(prior_dist, col=t_col(3,50) )
  
  
  #if(chain==2){
  #  
  #  legend("topright", legend=DEC_legend_names, fill =unlist(lapply(col_vec, function(c) t_col(c, 20))))
  #  
  #}
  
  
  
  polygon(density(DEC_post_dist[[chain]]), col=t_col("orange",20))
  
  
  # polygon(prior_dist, col=NA, border=t_col(3,50) )
  
  for( dir in 1:length(DEC_chains)){
    
    polygon(density(DEC_post_dist[[chain]]), col=NA, border = t_col("orange",0))
    
    
  }
} 


dev.off()   

#}


#}







###############plot

{
  
  pdf(paste("~/Projects/realfun_Biome/",dir_name,"/comp_density_rf_DEC_3bio.pdf", sep=""))
  
  
  
  par(mfrow=c(2,3))
  
  library(ggplot2)
  
  
  dat <- data.frame(Realized_Affinity_Loss_Rate = c(RFBS_post_dist[[1]]$X_rf2_l_s,post_dist[[2]]$X_l )
                    , lines = rep(c("RF real loss", "DEC loss"), each = length(post_dist[[1]]$X_rf2_l_s)))
  
  #ggplot(dat, aes(x = Realized_Affinity_Loss_Rate, fill = lines)) + geom_density(alpha = 0.5)
  
  ggplot(dat, aes(x = Realized_Affinity_Loss_Rate, fill = lines)) + geom_density(alpha = 0.5)+scale_x_continuous(limits = c(-0.1, 50)) +ggtitle("Realized Biome Affinity Losses")+  theme(plot.title = element_text(hjust = 0.5))   
  
  #plot(post_dist, main=chain_list[[chain]])
  
  #abline(v=sim_pars[sim,chain],col=2)
  
  dat <- data.frame(Realized_Affinity_Gain_Rate = c(post_dist[[1]]$X_rf2_g_s,post_dist[[2]]$X_g,post_dist[[1]]$X_rf2_g_d )
                    , lines = rep(c("RF real gain", "DEC gain", "RF real+fund gain"), each = length(post_dist[[1]]$X_rf2_g_s)))
  
  ggplot(dat, aes(x = Realized_Affinity_Gain_Rate, fill = lines)) + geom_density(alpha = 0.5)+scale_x_continuous(limits = c(-0.1, 20)) +ggtitle("Realized Biome Affinity Gains")+  theme(plot.title = element_text(hjust = 0.5))   
  
  
  
  #mtext(sim_pars_strings_unique[[sim]],                   # Add main title
  #      side = 3,
  #      line = - 2,
  #      outer = TRUE)
  #
  
  dev.off()
  
}










isin_HPD=sim_pars
ESS=sim_pars

post_mean=sim_pars
post_median=sim_pars





for (sim in 1:length(sim_pars_strings_unique)){
  #sort files based on same simulating pars
  runs=run_list[sim_pars_strings==sim_pars_strings_unique[[sim]]]
  
  test=lapply(runs, function(run) read.table(paste(dir_name,run, sep="/"), header = T))
  
  for (run in 1:length(runs)){
    
    chain_length=nrow(test[[run]])
    
    test[[run]] = test[[run]][(chain_length*burnin):chain_length,]
    
  }
  
  
  npars=(ncol(test[[1]])-4)/2
  
  chain_ind_vec=5:(4+npars)
  
  chain_list=colnames(test[[1]])[chain_ind_vec]
  
  par(mfrow=c(2,3))
  
  for (chain in 1:length(chain_ind_vec)){
    
    HPD=HPDinterval(as.mcmc(test[[run]][,chain_ind_vec[[chain]]]), prob=0.9)      
    ESS[sim,chain]=effectiveSize(as.mcmc(test[[run]][,chain_ind_vec[[chain]]]))   
    isin_HPD[sim,chain]=between(sim_pars[sim,chain],HPD[1],HPD[2])
    post_mean[sim,chain]=mean(test[[run]][,chain_ind_vec[[chain]]])
    post_median[sim,chain]=median(test[[run]][,chain_ind_vec[[chain]]])
    
    
    
  }
  
}



print(isin_HPD)

print(colSums(isin_HPD)/nrow(isin_HPD))
print(dir_name)



ESS_total=ESS
HPD_total=HPD
isin_HPD_total=isin_HPD



sum(ESS_total<100)/length(ESS_total)

post_mean_total=post_mean
post_median_total=post_median
sim_pars_total=sim_pars


{
  pdf(paste("~/Projects/realfun_Biome/",dir_name,"/posterior_median_lineplot.pdf", sep=""))
  
  par(mfrow=c(2,3))
  
  for (i in 1:ncol(post_median_total)){
    
    min=min(unlist(lapply(test, function(run) run[,chain_ind_vec[[chain]]]) ))
    min=min-min*.3
    max=max(unlist(lapply(test, function(run) run[,chain_ind_vec[[chain]]]) ))
    max=max+max*.3
    
    plot(c(0, 0), c(max, max), main=chain_list[[chain]],type = "n")
    
    lines(sim_pars_total[,i], post_median_total[,i])
    abline(a=0,b=1)
  }
  
  dev.off()
}

}
#  }

#}
#post_median_total/


#######

ESS_total=rbind(ESS_total,ESS)

HPD_total=rbind(HPD_total,HPD)
isin_HPD_total=rbind(isin_HPD_total,isin_HPD)
post_mean_total=rbind(post_mean_total,post_mean)
post_median_total=rbind(post_median_total,post_median)
sim_pars_total=rbind(sim_pars_total,sim_pars)
colSums(isin_HPD_total)/nrow(isin_HPD_total)

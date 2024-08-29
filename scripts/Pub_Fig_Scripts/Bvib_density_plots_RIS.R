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


new_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_3000000",list.files()  ))]
old_RFBS_dir_names=  list.files() [(grep("Bvib_3nB_Exp0p5_10000000",list.files()  ))][c( 5,1,6,7, 2,3,8,4)]

RFBS_dir_names=c(old_RFBS_dir_names, new_RFBS_dir_names)

evidence_vec=c("germination.only" ,"germination.conservative" ,"germination.bold", "leafing.conservative", "leafing.bold","USDA", "admat", "incf","excf"  )
name_vec=c("G_O" ,"G_C" ,"G_B", "L_C", "L_B","USDA", "A","I", "E" )

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

#RFBS_legend_order = c(9, 5,1,6,7, 2,3,8,4)

# Flatten the list into a character vector
RFBS_legend_names <- unlist(RFBS_legend_names)

new_RFBS_runs_df=cbind(RFBS_legend_names[9:length(RFBS_legend_names)], RFBS_dir_names_sort[9:length(RFBS_dir_names_sort)])

#new_RFBS_runs_df=cbind(RFBS_legend_names, RFBS_dir_names_sort)

# Sort by the number of dots (which corresponds to number of components after split) and then alphabetically

new_sorted_RFBS_jobs_df <- new_RFBS_runs_df[order(
  sapply(new_RFBS_runs_df[,1], function(x) length(strsplit(x, "\\.")[[1]])),
  new_RFBS_runs_df[,1]
),]


old_RFBS_runs_df=cbind(RFBS_legend_names[1:8], RFBS_dir_names_sort[1:8])

RFBS_legend_names

# Sort by the number of dots (which corresponds to number of components after split) and then alphabetically

old_sorted_RFBS_jobs_df <- old_RFBS_runs_df[order(
  sapply(old_RFBS_runs_df[,1], function(x) length(strsplit(x, "\\.")[[1]])),
  old_RFBS_runs_df[,1]
),]


old_RFBS_legend_names =match_string_elements_to_string_vec(string_vec    = old_sorted_RFBS_jobs_df[,1], 
                                                       string_pieces = c("E", "A", "I") , 
                                                       c("Old_Exc", "A", "Old_I"))



RFBS_legend_names=c(old_RFBS_legend_names, new_sorted_RFBS_jobs_df[,1])
RFBS_dir_names = c(old_sorted_RFBS_jobs_df[,2], new_sorted_RFBS_jobs_df[,2])

#RFBS_legend_names=c(new_sorted_RFBS_jobs_df[,1])
#RFBS_dir_names = c(new_sorted_RFBS_jobs_df[,2])


#RFBS_legend_names[[1]] = "Obs_only"
#RFBS_legend_names[[2]] = 
RFBS_legend_names[[length(RFBS_legend_names)+1]] = "Prior_only"


RFBS_Prior_dir="Bvib_3nB_Exp0p5_10000000_PO_admat_incf_excf_eco_allo_clado_2g_1g_1l_rf_gl_ds"

DEC_dir_names=  "Bvib_3nB_Exp0p5_1000000_DEC_eco_clado_2g_2l_gl"


setwd("/Volumes/michael.landis/Active/RFBS_RIS")

workdir=getwd()

#dir_name=dir_names[[1]]

tree=read.tree("viburnum_data_files/viburnum_sorted.tre")

rescale=max(nodeHeights(tree))
tree$edge.length

######make rfbs posterior objects

sim_dat=F
post_dist=list()
i=0

  
{
    
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
    DEC_file_list=lapply(DEC_dir_names, function(dir) list.files(dir))
    
    prior_file_list=list.files(RFBS_Prior_dir)
    
    
    #run_list=file_list[grep("__log.txt",unlist(file_list),fixed=FALSE)]
    #
    #
    prior_run=prior_file_list[grep("__log.txt",unlist(prior_file_list),fixed=FALSE)]
    #prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
    #
    RFBS_run_files=lapply(RFBS_file_list, function(file_list) file_list[grep("_log$",unlist(file_list),fixed=FALSE)])
   RFBS_prior_run=prior_file_list[grep("_log$",unlist(prior_file_list),fixed=FALSE)]
    RFBS_prior_test=read.table(paste(RFBS_Prior_dir,RFBS_prior_run[[1]], sep="/"), header = T)
    

    DEC_run_files=DEC_file_list[[1]][grep("_log$",unlist(DEC_file_list),fixed=FALSE)]
    #DEC_prior_run=prior_file_list[grep("_log$",unlist(prior_file_list),fixed=FALSE)]
    #DEC_prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
    

    
    burnin=0.5
####make rfbs post objects#####    
    
    RFBS_chains_full=lapply(1:length(RFBS_dir_names), function(dir) lapply(RFBS_run_files[[dir]], function(run) read.table(paste(RFBS_dir_names[[dir]],run, sep="/"), header = T) ))
    RFBS_chains=RFBS_chains_full
    npars=(ncol(RFBS_chains[[1]][[1]] )-4)/2
    chain_ind_vec=6:(5+npars)
    chain_list=colnames(RFBS_chains[[1]][[1]])[chain_ind_vec]
    RFBS_post_dist=list()
  
    {
    pdf(paste("EXP_0p5_bvib_RFBS_DEC_traces.pdf", sep=""),width = 8,height = 15)
    
      
for (dir in (1:length(RFBS_chains))[-106]){
  
  print(dir)
  
    RFBS_post_dist[[dir]]=list()
  
    for (run in 1:length(RFBS_chains[[dir]])){
      
      
      
      chain_length=nrow(RFBS_chains_full[[dir]][[run]])
      
      RFBS_chains[[dir]][[run]] = RFBS_chains_full[[dir]][[run]][(chain_length*burnin):chain_length,]
      

      
    }
      
    trace_cols <- distinctColorPalette(length(RFBS_chains[[dir]]))
    
    par(mfrow=c(3,2))
    
      
    for (chain in 1:length(chain_ind_vec)){
    
      
     # plot(RFBS_chains[[dir]][[run]][,chain_list[[chain]]]*rescale, type = "l", xlab = "", ylab="")
      for (run in 1:length(RFBS_chains[[dir]])){
    #    lines(RFBS_chains[[dir]][[run]][,chain_list[[chain]]], col=trace_cols[[run]])
      }
    

      RFBS_post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply( 1:length(RFBS_chains[[dir]]), function(run) RFBS_chains[[dir]][[run]][,chain_ind_vec[[chain]] ]))
    
      #lines(RFBS_prior_test[[chain_list[[chain]]]], col=t_col(3, 90))
      
   # post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply(test, function(run) run[,chain_ind_vec[[chain]] ]))
    
    }
    
    #title(  paste(RFBS_legend_names[[dir]]), line = -2, outer = TRUE)
    
      

     # post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply( RFBS_chains[[dir]], function(chain) run[,chain_ind_vec[[run]]]))
      
      
}
    
  dev.off()  
  
  
  
    }
    
#######make DEC posterior objects    
# 
#    DEC_chains_full=lapply( 1:length(DEC_run_files), function(run) read.table(paste(DEC_dir_names,DEC_run_files[[run]], sep="/"), header = T) )
#    DEC_chains=DEC_chains_full
#    DEC_npars=(ncol(DEC_chains[[1]] )-4)/2
#    DEC_chain_ind_vec=6:(5+DEC_npars)
#    DEC_chain_list=colnames(DEC_chains[[1]])[DEC_chain_ind_vec]
#    DEC_post_dist=list()
#    
#    
#    
#dir=1
#    #  DEC_post_dist[[dir]]=list()
#      
#      for (dir in 1:length(DEC_chains)){
#        
#        DEC_post_dist=list()
#        
#        for (run in 1:length(DEC_chains)){
#          
#          
#          
#          chain_length=nrow(DEC_chains_full[[run]])
#          
#          DEC_chains[[run]] = DEC_chains_full[[run]][(chain_length*burnin):chain_length,]
#          
#          
#          
#        }
#        
#        
#        for (chain in 1:length(DEC_chain_ind_vec)){
#          
#          
#          DEC_post_dist[[DEC_chain_list[[chain]]]]=unlist(lapply( 1:length(DEC_chains), function(run) DEC_chains[[run]][,DEC_chain_ind_vec[[chain]] ]))/rescale
#          
#          # post_dist[[chain_list[[chain]]]]=unlist(lapply(test, function(run) run[,chain_ind_vec[[chain]] ]))
#          
#        }
#        
#        
#        
#        # post_dist[[chain_list[[chain]]]]=unlist(lapply( DEC_chains, function(chain) run[,chain_ind_vec[[run]]]))
#        
#        
#      }
#      
#      
#      # post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply( DEC_chains[[dir]], function(chain) run[,chain_ind_vec[[run]]]))
#      
      
}
    

unique(as.numeric(RFBS_post_dist[[dir]][[chain]]))

RFBS_post_dens=list()

for( dir in (1:length(RFBS_post_dist))[-106]){
  RFBS_post_dens[[dir]]=list()
  for( chain in 1:length(RFBS_post_dist[[dir]])){
    RFBS_post_dens[[dir]][[chain]]=density(as.numeric(RFBS_post_dist[[dir]][[chain]]))
    
  }
  
names(RFBS_post_dens[[dir]])=   names(RFBS_post_dist[[dir]])

}

prior_chains=list()
    
######plot RFBS########
     
{

pdf(paste("3correctIEA_rescaled_EXP_0p5_bvib_RFBS_DEC_density.pdf", sep=""),width = 8,height = 10)


  
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
            "DEC and RFBS clado split prob"

            )

    for (z in 1:6){
      
      chain=c(2,1,5,3,4,6)[[z]]
      
      #plot(post_dist, main=chain_list[[chain]])
      prior_dist=density(unlist(sample(RFBS_prior_test[,chain_ind_vec[[chain]]]/rescale, nrow(RFBS_prior_test), replace = T)))
      
      #abline(v=sim_pars[sim,chain],col=2)
      x_min=min(unlist(lapply(RFBS_post_dist, function(i) min( i[[chain]]))))
      x_min=x_min-x_min*.3
      x_max=max(unlist(lapply(RFBS_post_dist, function(i) max( i[[chain]]))))
      if(chain!=length(chain_ind_vec)){
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
    




RFBS_post_dist=RFBS_post_dist[-106]  #exclude bad run
RFBS_legend_names=RFBS_legend_names[-106]  #exclude bad run
RFBS_chains = RFBS_chains[-106]


{
  
  
  
  pdf(paste("3biome_all_treatments_Donoghue_EXP_0p5_bvib_RFBS_density_violin.pdf", sep=""),width = 10,height = 10)
  
  
  
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
    par(mfrow = c(4, 1))
    
    par_names=c(
      "1->0",
      "0->1",
      "2->1",
      "0->2",
      "1->2",
      "clado split prob"
      
    )
    par_names=c(
      expression(italic(l)[1 %->% 0]), 
      expression(italic(g)[0 %->% 1]),
      expression(italic(l)[2 %->% 1]),
      expression(italic(g)[0 %->% 2]),
      expression(italic(g)[1 %->% 2]),
      expression(italic(b)[i]),
      expression(italic(w)[i]),
      expression(italic(e)[i])
    )
    
    
    for (z in 1:8){
      
      chain=c(2,1,5,3,4,6,7,8)[[z]]
      ymax= c(0.5, 2.5, 1.0, 1.5, 0.25, 1.0,1.0,1.0)
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

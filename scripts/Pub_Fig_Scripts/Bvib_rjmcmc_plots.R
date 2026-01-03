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

permutation_counts <- function(df, 
                               cols,
                               include_zero = FALSE,
                               keep_missing = TRUE) {
  stopifnot(all(cols %in% colnames(df)))
  library(data.table)
  
  ## 1.  Count every observed pattern ---------------------------------------
  DT      <- as.data.table(df[, cols, drop = FALSE])
  counts  <- DT[, .N, by = cols]                # N == frequency
  setorder(counts, -N)                          # biggest first
  
  ## 2.  Add unobserved permutations if requested ---------------------------
  if (keep_missing) {
    grid <- as.data.table(expand.grid(rep(list(c(0L, 1L)), length(cols)),
                                      KEEP.OUT.ATTRS = FALSE))
    setnames(grid, cols)
    counts <- merge(grid, counts, by = cols, all.x = TRUE)
    counts[is.na(N), N := 0L]
  }
  
  ## 3.  Optionally drop the all-zero row -----------------------------------
  if (!include_zero) {
    counts <- counts[rowSums(counts[, ..cols]) > 0]
  }
  
  counts[]
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
  #"Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss_RJa",
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


evidence_vec=c( "Bvib_3nB_Exp0p5_1000000_PO_Foss_noincf_cons_noexcf_3.biomes.germination.only.germination.bold.leafing.bold_eco_allo_clado_2g_1g_1l_2sw_rf_gl_ds_Foss",
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

RFBS_legend_names[RFBS_legend_names == ""] = "None"

#RFBS_legend_order = c(9, 5,1,6,7, 2,3,8,4)

# Flatten the list into a character vector
RFBS_legend_names <- unlist(RFBS_legend_names)

RFBS_legend_names  = gsub(".", "_", RFBS_legend_names, fixed = T )

new_RFBS_runs_df = cbind(RFBS_legend_names, RFBS_dir_names_sort)

new_RFBS_runs_df = new_RFBS_runs_df[!duplicated(new_RFBS_runs_df[,1]), ]


RFBS_treatments_sorted=c("Prior", "None", "A", "Ic","Ib", "Ec", "Eb", "A_Ic", "A_Ib", "A_Ec", "A_Eb","Ic_Ec", "Ic_Eb", "Ib_Ec", "Ib_Eb", "A_Ic_Ec", "A_Ic_Eb", "A_Ib_Ec", "A_Ib_Eb")

sorted_RFBS_runs_df= new_RFBS_runs_df[match(RFBS_treatments_sorted,new_RFBS_runs_df[,1] ) , ]

RFBS_dir_names = sorted_RFBS_runs_df[-1,2]
RFBS_Prior_dir = sorted_RFBS_runs_df[1,2]

RFBS_legend_names = sorted_RFBS_runs_df[-1,1]


setwd("/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp/viburnum/resub")

workdir="/Volumes/michael.landis/Active/Sean/RFBS/outfiles/emp/viburnum/resub/"

RFBS_dir_names = paste0(workdir, RFBS_dir_names)
RFBS_Prior_dir = paste0(workdir, RFBS_Prior_dir)

#dir_name=dir_names[[1]]

tree=read.tree("/Volumes/michael.landis/Active/Sean/RFBS/data/emp/viburnum_data_files/viburnum_sorted.tre")

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

prior_file_list=list.files(RFBS_Prior_dir)
#
prior_run=prior_file_list[grep("__log.txt",unlist(prior_file_list),fixed=FALSE)]
#prior_test=read.table(paste(prior_dir,prior_run, sep="/"), header = T)
#
RFBS_run_files=lapply(RFBS_file_list, function(file_list) file_list[grep("_log$",unlist(file_list),fixed=FALSE)])
RFBS_prior_run=prior_file_list[grep("_log$",unlist(prior_file_list),fixed=FALSE)]
RFBS_prior_test=read.table(paste(RFBS_Prior_dir,RFBS_prior_run[[1]], sep="/"), header = T)

burnin=0.5
####make rfbs post objects#####    
######this loads the chains takes awhile and need connection to RIS
RFBS_chains_full=lapply(1:length(RFBS_dir_names), function(dir) lapply(RFBS_run_files[[dir]], function(run) read.table(paste(RFBS_dir_names[[dir]],run, sep="/"), header = T) ))
RFBS_chains_prior= lapply(RFBS_run_files[[1]], function(run) read.table(paste(RFBS_dir_names[[1]],run, sep="/"), header = T) )

###############################################


{
  
  RFBS_chains=RFBS_chains_full
  npars=19 #(ncol(RFBS_chains[[1]][[1]] )-4)/2
  
  
  #chain_ind_vec=6:(5+npars)
  
  #chain_ind_vec = c(6:11, 18:20)
  chain_ind_vec = c(6:20)
  
  chain_list=colnames(RFBS_chains[[1]][[1]])[chain_ind_vec]
  RFBS_post_dist=list()
  
  {
    for (dir in (1:length(RFBS_chains))){
      print(dir)
      RFBS_post_dist[[dir]]=list()
      for (run in 1:length(RFBS_chains[[dir]])){
        chain_length=nrow(RFBS_chains_full[[dir]][[run]])
        #RFBS_chains[[dir]][[run]] = RFBS_chains_full[[dir]][[run]][(chain_length*burnin):chain_length,]
        RFBS_chains[[dir]][[run]] = RFBS_chains_full[[dir]][[run]][5,]
        
      }
      for (chain in 1:length(chain_ind_vec)){
        RFBS_post_dist[[dir]][[chain_list[[chain]]]]=unlist(lapply( 1:length(RFBS_chains[[dir]]), function(run) RFBS_chains[[dir]][[run]][,chain_ind_vec[[chain]] ]))
      }
    }
  }
}


chain_par_vec = c(6:11, 18:20)-5
chain_par_names = chain_list[chain_par_vec]
chain_rj_vec = c(12:17) - 5

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
    RFBS_post_dens[[dir]][[i]]      = density(as.numeric(RFBS_post_dist[[dir]][[RFBS_names[chain]]]))
    names(RFBS_post_dens[[dir]][i]) = names(RFBS_post_dist[[dir]])[chain]
    
  }
  #names(RFBS_post_dens[[dir]])=   names(RFBS_post_dist[[dir]])[chain_par_vec]
  for(  i in 1:length(chain_rj_vec)){
    chain = chain_rj_vec[[i]]
    
    RFBS_rj_dens[[dir]][[i]]=density(as.numeric(RFBS_post_dist[[dir]][[RFBS_names[chain]]]))
    names(RFBS_rj_dens[[dir]][i])=   names(RFBS_post_dist[[dir]])[chain]
    
  }
}

prior_chains=list()

chain_rj_names = c("switch_rf1_l_s",  
                   "switch_rf2_l_d",   "switch_rf1_g_s",  
                   "switch_rf12_sw_s", "switch_rf2_l_s",  
                   "switch_rf2_g_d",   "switch_rf2_g_s"  )
rj_post_array_list = lapply(c(1, 15, 18), function(dir) do.call(cbind,RFBS_post_dist[[dir]][chain_rj_names]))
names(rj_post_array_list)  = unlist(lapply(c(2, 16, 19), function(dir) RFBS_legend_names[[dir]]))



library(tools)
library(coda)
library(ggplot2)
library(patchwork)
library(RevGadgets)
library(RColorBrewer)

# Directory Setup ---------------------------------------------------------


# Basic Jointplot ---------------------------------------------------------

#data <- read.csv(file=paste(data_dir,"/concatenated.model.log",sep=""),sep="\t",header=TRUE)

for (experiment  in 1:length(rj_post_array_list)){
{
  data = rj_post_array_list[[experiment]]
  
  experi_name = names(rj_post_array_list)[[experiment]]
  
  
  colnames(data)
  
  joint <- function(param1, param2) {
    param1b <- param1#paste("rj_",gsub("\\[|\\]", ".", param1),sep="")
    param2b <- param2 #paste("rj_",gsub("\\[|\\]", ".", param2),sep="")
    rj1 <- data[,param1b]
    rj2 <- data[,param2b]
    joints <- c(0,0,0,0)
    for (i in 1:length(rj1)) {
      if (rj1[i] == 1) {
        if (rj2[i] == 1) {joints[1] <- joints[1] + 1}
        else {joints[2] <- joints[2] + 1}
      }
      else {
        if (rj2[i] == 1) {joints[3] <- joints[3] + 1}
        else {joints[4] <- joints[4] + 1}
      }
    }
    joints <- round(joints/sum(joints),2)
    return(joints)
  }
  
  jointplot <- function(param1,param2, xname, yname, title) {
    joints <- joint(param1,param2)
    dataframe <- data.frame(matrix(c(0,1,2,0,1,2),ncol=2))
    plot <- ggplot(dataframe,aes(x=X1,y=X2)) +
      geom_rect(xmin=1,xmax=2,ymin=1,ymax=2,color="black",fill="blue",alpha=joints[1]) +
      geom_rect(xmin=1,xmax=2,ymin=0,ymax=1,color="black",fill="blue",alpha=joints[2]) +
      geom_rect(xmin=0,xmax=1,ymin=1,ymax=2,color="black",fill="blue",alpha=joints[3]) +
      geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,color="black",fill="blue",alpha=joints[4]) +
      #annotate(geom="text",x=1.5,y=1.5,label=joints[1], size=10, colour = "white") +
      #annotate(geom="text",x=1.5,y=.5,label=joints[2] , size=10, colour = "white") +
      #annotate(geom="text",x=.5,y=1.5,label=joints[3] , size=10, colour = "white") +
      #annotate(geom="text",x=.5,y=.5,label=joints[4]  , size=10, colour = "white") +
      annotate(geom="text",x=1.5,y=1.5,label=joints[1], size=10) +
      annotate(geom="text",x=1.5,y=.5,label=joints[2], size=10) +
      annotate(geom="text",x=.5,y=1.5,label=joints[3], size=10) +
      annotate(geom="text",x=.5,y=.5,label=joints[4], size=10) +
      
      scale_x_continuous(limits=c(0,2),breaks=c(.5,1.5),labels=c("Off","On")) +
      scale_y_continuous(limits=c(0,2),breaks=c(.5,1.5),labels=c("Off","On")) +
      labs(title=title,x=xname,y=yname) +
      theme_classic() +
      theme(aspect.ratio=1,plot.title=element_text(hjust=.5),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(angle=90,hjust=0.5,vjust= -.3),axis.text.x=element_text(vjust=1.5))
    return(plot)
  }
  
  
  par_names=c(
    expression(paste("Enabled loss "      , italic(l)[1 %->% 0])), 
    expression(paste("Enabled gain "      , italic(g)[0 %->% 1])),
    expression(paste("Established switch ", italic(sw)[2])),
    expression(paste("Established loss "  , italic(l)[2 %->% 1])),
    expression(paste("Double gain "       , italic(g)[0 %->% 2])),
    expression(paste("Established gain "  , italic(g)[1 %->% 2])),
    expression(paste("Double loss "       , italic(l)[2 %->% 0]))
    
    #expression(paste("Speciation event "      , italic(b))          ),
    #expression(paste("Speciation event "      , italic(s))          ),
    #expression(paste("Speciation event "      , italic(e))          )
  )
  {
  w1 <- jointplot("switch_rf1_g_s","switch_rf2_g_d", par_names[[2]], par_names[[5]], "Enabled Gain vs.\nDouble Gain")
  w2 <- jointplot("switch_rf2_g_s","switch_rf2_g_d", par_names[[6]], par_names[[5]], "Established Gain vs.\nDouble Gain")
  w3 <- jointplot("switch_rf1_g_s","switch_rf2_g_s", par_names[[2]], par_names[[6]], "Enabled Gain vs.\nEstablished Gain")
  w4 <- jointplot("switch_rf1_g_s","switch_rf12_sw_s", par_names[[2]], par_names[[3]], "Enabled Gain vs.\nEstablished Switch")
  w5 <- jointplot("switch_rf2_g_s","switch_rf12_sw_s", par_names[[6]], par_names[[3]], "Established Gain vs.\nEstablished Switch")
  w6 <- jointplot("switch_rf2_g_d","switch_rf12_sw_s", par_names[[5]], par_names[[3]], "Double Gain vs.\nEstablished Switch")

  
  
  
  axis_plot <- ggplot() +
    #labs(x=bquote("Categorical (" ~ sigma ~ ")"),y=bquote("Quantitative (" ~ phi ~ ")")) +
    theme_classic() +
    theme(aspect.ratio=1.5,line=element_blank())
  inner_plot <- w1 + w2 + w3 + w4 + w5 + w6 +
    plot_layout(ncol=2) &
    theme(text=element_text(size=))
  joint_plot <- axis_plot + inset_element(inner_plot,left=0,bottom=0,right=1,top=1)
  
  pdf(file=paste("~/Projects/RFBS-main/outfiles/emp/viburnum/resub_figs/gain_RJ_pairjoint_2l_", experi_name,  ".pdf", sep=""), width = 7.5, height =  10)
  print(joint_plot)
  dev.off()
  }
  
  {
  w1 <- jointplot("switch_rf1_l_s","switch_rf2_l_d"  , par_names[[1]], par_names[[7]], "Enabled Loss vs.\nDouble Loss")
  w2 <- jointplot("switch_rf2_l_s","switch_rf2_l_d"  , par_names[[4]], par_names[[7]], "Established Loss vs.\nDouble Loss")
  w3 <- jointplot("switch_rf1_l_s","switch_rf2_l_s"  , par_names[[1]], par_names[[4]], "Enabled Loss vs.\nEstablished Loss")
  w4 <- jointplot("switch_rf1_l_s","switch_rf12_sw_s", par_names[[1]], par_names[[3]], "Enabled Loss vs.\nEstablished Switch")
  w5 <- jointplot("switch_rf2_l_s","switch_rf12_sw_s", par_names[[4]], par_names[[3]], "Established Loss vs.\nEstablished Switch")
  w6 <- jointplot("switch_rf2_l_d","switch_rf12_sw_s", par_names[[7]], par_names[[3]], "Established Loss vs.\nEstablished Switch")
  
  
  
  
  axis_plot <- ggplot() +
    #labs(x=bquote("Categorical (" ~ sigma ~ ")"),y=bquote("Quantitative (" ~ phi ~ ")")) +
    theme_classic() +
    theme(aspect.ratio=1.5,line=element_blank())
  inner_plot <- w1 + w2 + w3 + w4 + w5 + w6 +
    plot_layout(ncol=2) &
    theme(text=element_text(size=))
  joint_plot <- axis_plot + inset_element(inner_plot,left=0,bottom=0,right=1,top=1)
  
  pdf(file=paste("~/Projects/RFBS-main/outfiles/emp/viburnum/resub_figs/loss_RJ_pairjoint_2l_", experi_name,  ".pdf", sep=""), width = 7.5, height =  10)
  print(joint_plot)
  dev.off()
  }
}
}
  # Big Jointplot -----------------------------------------------------------
  for (experiment  in 1:length(rj_post_array_list)){
    
  
  par_names=c(
    expression(paste("Enabled loss "      , italic(l)[1 %->% 0])), 
    expression(paste("Double loss "       , italic(l)[2 %->% 0])),
    expression(paste("Enabled gain "      , italic(g)[0 %->% 1])),
    expression(paste("Established switch ", italic(sw)[2])),
    expression(paste("Established loss "  , italic(l)[2 %->% 1])),
    expression(paste("Double gain "       , italic(g)[0 %->% 2])),
    expression(paste("Established gain "  , italic(g)[1 %->% 2]))

    #expression(paste("Speciation event "      , italic(b))          ),
    #expression(paste("Speciation event "      , italic(s))          ),
    #expression(paste("Speciation event "      , italic(e))          )
  )
  
  data_onoff <- data
  rjparams <- colnames(data)
  
  sums <- rowSums(data_onoff)
  onoff_sums <- cbind(data_onoff, sums)
  sums <- rep(sums, 10)
  onoff <- data.frame(sums)
  total <- nrow(data_onoff)
  
  total_df <- data.frame(matrix(ncol=2,nrow=1))
  for (i in rjparams) {
    param <- i
    count <- sum(data_onoff[,param])
    row <- data.frame(t(c(i,count)))
    total_df <- rbind(total_df,row)
  }
  total_df <- total_df[-(1),]
  colnames(total_df) <- c("param","percent")
  total_df$percent <- as.numeric(total_df$percent)/total
  
  stacked_df <- data.frame(matrix(ncol=3,nrow=1))
  for (i in 0:6) {
    onoff_sums_i <- onoff_sums[which(onoff_sums[,ncol(onoff_sums)]==i),]
    for (j in rjparams) {
      param <- j
      count <- sum(onoff_sums_i[,param])
      row <- data.frame(t(c(i,j,count)))
      stacked_df <- rbind(stacked_df,row)
    }
  }
  stacked_df <- stacked_df[-(1),]
  colnames(stacked_df) <- c("onoff","param","count")
  stacked_df$count <- as.numeric(stacked_df$count)
  
  levels <- colnames(data)
  labels <- par_names
  
  
  palette1 <- brewer.pal(8,"Dark2")
  palette2 <- brewer.pal(8,"Set2")
  palette3 <- c(rbind(palette1))
  
  total_plot <- ggplot(total_df,aes(x=factor(param,levels=levels),y=percent,fill=factor(param,levels=levels))) +
    geom_bar(stat="identity") +
    scale_fill_manual(labels=labels,values=palette3) +
    scale_y_continuous(limits=c(0,1),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0)) +
    labs(fill="Parameter",x="Parameter",y="Frequency") +
    theme_bw() +
    theme(aspect.ratio=.4,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(8,8,4,8))
  
  onoff_plot <- ggplot(onoff, aes(sums)) +
    geom_histogram(aes(y=..density..), binwidth=1, boundary=.5, color="black", size=.1, fill="white") +
    geom_density(aes(linetype="Observed"), adjust=4, linewidth=1) +
    geom_density(aes(rbinom(nrow(onoff),6,.5), linetype="Prior"), adjust=4, linewidth=1) +
    scale_linetype_manual(values=c("Observed"="solid","Prior"="dashed")) +
    scale_x_continuous(breaks=seq(0,6,1), limits=c(-.5,6.5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,2),breaks=seq(0.0, 1.5, .5),labels=paste(seq(0.0, 1.5, .5)),expand=c(0,0)) +
    labs(y="Density", x=NULL, linetype="Distribution") +
    theme_bw() +
    theme(aspect.ratio=.3,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.margin=margin(4,8,4,8))
  
  stacked_plot <- ggplot(stacked_df,aes(x=onoff,y=count,fill=factor(param,levels=levels))) +
    geom_bar(stat="identity",position="fill",width=1,color="black",size=.1) +
    scale_x_discrete(limits=factor(seq(0,6,1)),expand=c(0,0)) +
    scale_y_continuous(breaks=seq(0,1,.25),labels=c("0.00","0.25","0.50","0.75","1.00"),expand=c(0,0)) +
    scale_fill_manual(labels=labels,values=palette3) +
    labs(fill="Parameter",x="Number of 'ON' Parameters",y="Representation By Bin") +
    theme_bw() +
    theme(aspect.ratio=.75,panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none",plot.margin=margin(4,8,8,8))
  
  big_onoff_plot <- total_plot + onoff_plot + stacked_plot + plot_layout(ncol=1,guides="collect")
  
  pdf(file=paste("~/Projects/RFBS-main/outfiles/emp/viburnum/resub_figs/2l_big_onoff_plot_",experi_name,".pdf", sep=""))
  print(big_onoff_plot)
  dev.off()
  
}
  




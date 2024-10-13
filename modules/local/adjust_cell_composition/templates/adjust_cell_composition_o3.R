#!/usr/bin/env Rscript
##This script reads M-values and Beta values of methylation arrays and corrects for cell composition (samples from whole blood)
##Use after visualizing cell composition in the different samples and deciding if you want to correct for that effect
##Input = Beta values
##Output = M-values and beta values corrected for cell commpositions

library(readr)
library(ChAMP)
library(lumi)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
bVals <- as.matrix(read_csv("$bVALS_SNPPROBES") %>% tibble::column_to_rownames("probe"))

##Correct for cell composition using ChAMP 
bVals_corrected <- champ.refbase(beta = bVals, arraytype = "EPIC")
bVals <- bVals_corrected\$CorrectedBeta
mVals <- beta2m(bVals)

##Save all necessary R objects for later use
write_csv(as.data.frame(mVals), "cmVals.cell_comp.csv")
write_csv(as.data.frame(bVals), "cbVals.cell_comp.csv")
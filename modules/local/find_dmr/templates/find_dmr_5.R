#!/usr/bin/env Rscript
##This script reads Beta values of methylation arrays and computes DMRs, Differentially methylated regionss
##Input = beta values of pre-processed methylation data 
##Input = metadata information with sample ids and group information (categories: case, control or different groups)
##output = csv files with p values of DMRs found by the chosen method

library(dplyr)
library("readr")
library("minfi")
library("ChAMP")
library(DMRcate)
library(rio)
library(readr)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd, cell_composition_correction.R and rem_conf_probes_adj_age.R
###Load necessary data
bVals <- read_csv("$bVALS_SNPPROBES")
metadata <- read_csv("$extensive_metadata") %>%
  filter(sample_id %in% c(colnames(bVals)))
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(bVals)), metadata\$sample_id),]

##Choose the samples field
##Choosing only the field with information of the category of each samples (can be multiple)
Samples <- metadata\$sample_id
Class <- metadata\$group

##Choose the method: can be "Bumphunter", "ProbeLasso" or "DMRcate" but it actually doesn't work for DMRcate (not maintained anymore through ChAMP, 
###so I calculate it with its explicit function here)
Method = "Bumphunter"

##Choose array type: ("EPIC" or "450K")
ARRAY = "EPIC"

###Choose adjusted P value
P = 0.05
P = 1 # NOTE: for development

###Choose number of minimum probes by region
Min = 5

###Choose maximum gap between 2 probes in a region
Max = 300 
Max = 10000 # NOTE: for development

### Choose number of bootstraps
B = 1000

###Choose bootstrapping method can be ‘bootstrap’ or ‘permutation’
Boot = "bootstrap"

##Define the number of cores
N = 1

###Defininf the DMRcate function
DMRcate_manual_run = function(beta_matrix,
                              pheno,
                              arraytype,
                              dist = 2, 
                              mafcut = 0.05,
                              fdr = 1,
                              lambda = 300,
                              C = 2)
{
  require(DMRcate)
  myMs <- logit2(beta_matrix)
  myMs <- rmSNPandCH(myMs, dist = dist, mafcut = mafcut)
  design <- model.matrix(~ pheno)
  if(arraytype == "EPIC"){
    myannotation <- cpg.annotate(datatype = "array",
                                 fdr = 1, 
                                 myMs, 
                                 design = design,
                                 coef = ncol(design), 
                                 analysis.type = "differential",
                                 annotation = c(array = "IlluminaHumanMethylation450k", 
                                                annotation = "ilmn12.hg19"), 
                                 what = "M")
  } else {
    myannotation <- cpg.annotate(datatype = "array",
                                 fdr = 1, 
                                 myMs, 
                                 design = design,
                                 coef = ncol(design), 
                                 analysis.type = "differential",
                                 annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                                annotation = "ilm10b4.hg19"), 
                                 what = "M")
  }
  dmrcoutput <- dmrcate(myannotation, lambda = lambda, C = C, ) 
  DMR <- as.data.frame(extractRanges(dmrcoutput, genome = "hg19"))
  rownames(DMR) <- paste("DMR", 1:nrow(DMR), sep="_")
  message("DMRcate detected ", nrow(DMR)," DMRs with mafcut as = ", mafcut, ".")
  if(nrow(DMR) == 0) 
    stop("No DMR detected.")
  
  return(DMR)
}

####Running a method depending on the user's choice

if (Method %in% c("Bumphunter", "ProbeLasso")) {
  
  ###Compute Champ's bumphunter based DMR permutations
  dmr_champ <- champ.DMR(beta = as.matrix(bVals[,Samples]),
                         pheno = Class,
                         arraytype = ARRAY,
                         method = Method,
                         minProbes = Min,
                         adjPvalDmr = P,
                         maxGap = Max,
                         B = B,
                         nullMethod = Boot,
                         cores = N) 
  
  export_list(dmr_champ, file = "dmr_champ.%s.csv")
}

###Computing DMRcate with the explicit formula, as it doesn't work anymore in ChAMP
if (Method == "DMRcate") {
  ##DMRcate computations
  dmrcate <- DMRcate_manual_run(bVals[,Samples], 
                                Class, 
                                ARRAY)
  write_csv(dmrcate, "dmr_dmrcate.csv")
}
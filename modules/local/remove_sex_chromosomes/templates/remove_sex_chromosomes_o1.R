#!/usr/bin/env Rscript

##This script does the preprocessing of methylation data  to normalize and filter the data including sample QC, probes QC, removing sex chromosomes, SNPs probes, cross-reactive
#probes and correcting for cell compostition. No batch effect have been noted to require
#correction (see pre-processing.Rmd). Finally, the script calculates methylation values.
##input= Raw data .idat files
##output=normalized  M and beta values (.RData format)

##Requirements
library("readr")
library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(tibble)
library(FlowSorted.Blood.EPIC)

##Constants
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

##Extract list of cross-reactive probes found using the script xreactive_probes_finder.R
#xprobes <- read_csv("x_reactive_probes.csv") # is this needed?

##Load data (for more detailed descriptions, see pre-processing.Rmd)
get(load("$RData_SNPPROBES"))
get(load("$RData_rgSet"))
#"mSetSqFlt.RData"
#"rgSet.RData"


##Remove probes from X and Y chromosomes
keep <- !(featureNames(mSetSqFlt) %in% annotation\$Name[annotation\$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]

##Remove probes with known SNPs, default is to remove all, but you can define a threshold of highest frequency of the minor allele
snps <- getSnpInfo(mSetSqFlt) #GRset
mSetSqFlt <- addSnpInfo(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps=c("SBE","CpG", "Probe"), maf=0) ##Can add "Probe" to the SNP list and can change maf threshold



##Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

# NOTE: Why is this here?
##Calculate cell composition estimation
cellCounts <- estimateCellCounts(rgSet) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_id") %>% 
  write_csv("cell_count.csv")

##Save all necessary R objects for later use
save(mVals, file = "mVals.RData")
save(bVals, file = "bVals.RData")
save(mSetSqFlt, file = "mSetSqFlt.RData")
save(rgSet, file = "rgSet.RData") # NOTE: Not needed?

library(dplyr)
write_csv(as.data.frame(bVals), "mVals.csv")
write_csv(as.data.frame(mVals), "bVals.csv")
write_csv(as.data.frame(bVals) %>% tibble::rownames_to_column(var = 'probe') , "mVals.csv")
write_csv(as.data.frame(mVals) %>% tibble::rownames_to_column(var = 'probe'), "bVals.csv")

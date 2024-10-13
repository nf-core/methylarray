#!/usr/bin/env Rscript
##This script reads Beta values of methylation arrays and computes blocks of dna regions with many DMPs, 
##Input = beta values of pre-processed methylation data 
##Input = metadata information with sample ids and group information (categories: case, control or different groups)
##output = csv files with p values of blocks 

library(dplyr)
library("readr")
library("ChAMP")
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

##Choose array type: ("EPIC" or "450K")
ARRAY = "EPIC"

###Choose number of minimum number of regions per block
Min = 5

###Choose maximum gap between 2 clusters in a block
Max = 250000

###Computing blocks
blocks <- champ.Block(beta = bVals[,Samples], 
                       pheno = Class, 
                       arraytype = ARRAY, 
                       maxClusterGap = Max, 
                       minNum = Min)

  write_csv(blocks, "blocks_champ.csv")
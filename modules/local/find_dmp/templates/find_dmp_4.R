#!/usr/bin/env Rscript
##This script reads Beta values of methylation arrays and computes DMPs, Differentially methylated positions
##Input = beta values of pre-processed methylation data 
##Input = metadata information with sample ids and group information (categories: case, control or different groups)
##output = csv files with p-values for the DMPs in each comparison of categories

##Computations of DMPs using ChAMP and minfi

library(dplyr)
library("readr")
library("minfi")
library(tibble)
library("ChAMP")
library(rio)
library(readr)


##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd, cell_composition_correction.R and rem_conf_probes_adj_age.R
bVals <- read_csv("$bVALS_SNPPROBES")
metadata <- read_csv("$extensive_metadata") %>%
            filter(sample_id %in% c(colnames(bVals)))
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(bVals)), metadata\$sample_id),]

##Choose the samples field
##Choosing only the field with information of the category of each samples (can be multiple)
Samples <- metadata\$sample_id
Class <- metadata\$group

##Choose the adjustment method: can be "BH", 
Method = "BH"

##Choose array type: ("EPIC" or "450K")
ARRAY = "EPIC"

###Choose adjusted P value
P = 0.05

###Computing DMPs
dmp_data <- champ.DMP(beta = bVals[,Samples], 
          pheno = Class, 
          adjPVal = P, 
          adjust.method = Method, 
          arraytype = ARRAY) 

export_list(dmp_data, file = "dmp_champ.%s.csv")

###Finding DMPs with another method (required binary categories of "Class", it will work with multiple categories but the do not specify which comparison it is)
dmp_minfi <- dmpFinder(bVals[,Samples], 
                 Class, 
          type = "categorical") %>% 
          filter(pval < P) %>%
          arrange(pval) 

write_csv(dmp_minfi, file = "dmp_minfi.csv")
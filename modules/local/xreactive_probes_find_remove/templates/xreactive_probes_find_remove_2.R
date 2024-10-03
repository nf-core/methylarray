###!usr/R
### This script finds cross reactive probes in IlluminaEPIC array and removes them, It is also valid for all types of Illumina arrays
### input = normalized  and filtered mSet object 
##output = mSet object (.RData format), normalized  M and beta values (csv format)

library("minfi")
library("IlluminaHumanMethylationEPICmanifest")
library(limma)
library(missMethyl)
library("IlluminaHumanMethylation450kmanifest")
library(minfiData)
library(DMRcate)
library(DNAmCrosshyb)
library(argparse)
library("dplyr")
library("readr")

# Create argument parser
parser <- ArgumentParser(description = 'Preprocess methylation data')
parser$add_argument('rdata', type = 'character', help = '.RData from preprocess')
parser$add_argument('bs_genome_path', type = 'character', help = 'bs_genome_path')

# Parse the command-line arguments
args <- parser$parse_args()

# Get the input arguments
RData_PREPROCESSING <- args$rdata
Genome_PATH <- args$bs_genome_path

####The human genome version should be downloaded to data/genome_bs and here we specify the file 
#PATH <- "data/genome_bs/hg19"
PATH <- Genome_PATH

###Choose minimum width, maximum width and the step for the probes cross reaction
MIN = 30
MAX = 50
STEP = 5

##Read prestored data (see pre-processing.Rmd for more details)
get(load(RData_PREPROCESSING))

##Cross-reactive probes detection, 
probes <- rownames(mSetSqFlt)
matches <- map_probes(probes, path = PATH, 
                      # chromosomes was 'all' but is reduced to 2 for testing
                      array = "EPIC", chromosomes = "2", min_width = MIN, max_width = MAX, step_size = STEP, 
                      allow_mismatch = FALSE, allow_INDEL = FALSE, verbose = TRUE)

###Identify probes mapping to multiple sites on the genome
nr_matches <- get_nr_matches_per_probe(matches)
nr_matches %>% 
  data.frame() %>% 
  ###Choose minimum number of bp for the mismatch to be retained: 30, 35, 40, 45 or 50
  filter(bp30 != "1" | bp35 != "1" | bp40 != "1" | bp45 != "1" | bp50 != "1") %>% 
  select(Probe) %>%
  write_csv("x_reactive_probes.csv")

rep_matches <- find_repeat_overlaps(matches, genome_build = "hg19", min_overlap = "any") 
  rep_matches %>%      ####manually check the number of probes to decide he width threshold you want to have
    filter(repeat_overlap == "TRUE") %>% 
    filter(width == MAX) %>% 
    distinct(Probe)
  
nr_matches <- as.data.frame(nr_matches)
  
##Remove cross-reactive probes from list I made
keep <- !(featureNames(mSetSqFlt) %in% nr_matches$Probe)
mSetSqFlt <- mSetSqFlt[keep,] 

##Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

##Save all necessary R objects for later use
write_csv(as.data.frame(mVals), "mVals_noXprob.csv")
write_csv(as.data.frame(bVals), "bVals_noXprob.csv")
save(mSetSqFlt, file = "mSetSqFlt_noXprob.RData")



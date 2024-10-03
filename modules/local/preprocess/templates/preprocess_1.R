#!/usr/bin/env Rscript

# Load necessary libraries
library(argparse)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)

# Create argument parser
parser <- ArgumentParser(description = 'Preprocess methylation data')
parser$add_argument('idat_folders', type = 'character', help = 'Path to the folder containing .idat files')
parser$add_argument('samplesheet_name', type = 'character', help = 'Name of the sample sheet (CSV format)')

# Parse the command-line arguments
args <- parser$parse_args()

# Get the input arguments
dataDirectory <- args$idat_folders
metharray_sheet <- args$samplesheet_name

# Set constants
P = 0.01
Norm_method = "preprocessQuantile"  # Choose normalization method: "preprocessFunnorm" or "preprocessQuantile"
Sample_Name <- "Sample_Name"        # Column name that encodes the sample names

# Print some info for debugging
cat("Data directory:", dataDirectory, "\n")
cat("Sample sheet:", metharray_sheet, "\n")

# Load data
list.files(dataDirectory, recursive = TRUE)

# Reading the sample sheet in CSV format
targets <- read.metharray.sheet(dataDirectory, pattern = metharray_sheet)
rgSet <- read.metharray.exp(targets = targets)

# Assign significant names to the columns (samples)
targets$ID <- targets[[Sample_Name]]
sampleNames(rgSet) <- targets$ID

# QC: Detection p-values of the signal quality
detP <- detectionP(rgSet)
keep <- colMeans(detP) < P
rgSet <- rgSet[, keep]

# Remove poor quality samples from targets and detection p-value table
targets <- targets[keep, ]
detP <- detP[, keep]

# Normalize the data
if (Norm_method == "preprocessQuantile") {
  mSetSq <- preprocessQuantile(rgSet)
} else if (Norm_method == "preprocessFunnorm") {
  mSetSq <- preprocessFunnorm(rgSet)
}

# Filtering: Put probes in the same order in mSetSq and detP
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]

# Filter probes failing in 1 or more samples (change threshold if needed)
keep <- rowSums(detP < P) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep, ]

# Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

# Save results
write.csv(as.data.frame(mVals), "mVals.csv")
write.csv(as.data.frame(bVals), "bVals.csv")
save(mSetSqFlt, file = "mSetSqFlt.RData")
save(rgSet, file = "rgSet.RData")

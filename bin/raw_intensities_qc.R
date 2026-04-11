#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rgset", type="character", default=NULL, help="Path to RGset RDS file")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path <- opt$rgset
if (is.null(rg_path) || !file.exists(rg_path)) stop("RGset RDS not found: ", rg_path)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

req("minfi")

tryCatch({
  say("== Raw intensities QC (minfi) ==")
  rg0 <- readRDS(rg_path)

  ms_raw <- minfi::preprocessRaw(rg0)
  qc <- minfi::getQC(ms_raw)

  qc_df <- data.frame(Sample = rownames(qc), mMed = qc$mMed, uMed = qc$uMed, stringsAsFactors = FALSE)
  write.table(qc_df, file = "qc_raw.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  low <- rownames(qc)[qc$mMed < 10.5 | qc$uMed < 10.5]
  low_df <- data.frame(
    Sample = low,
    Meth_median = qc$mMed[low],
    Unmeth_median = qc$uMed[low],
    stringsAsFactors = FALSE
  )
  write.table(low_df, file = "low_intensity_samples.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  save_png("01_plotQC_raw_intensities.png", minfi::plotQC(qc, badSampleCutoff = 10.5))

  say("Wrote: 01_plotQC_raw_intensities.png, qc_raw.tsv, low_intensity_samples.tsv")
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

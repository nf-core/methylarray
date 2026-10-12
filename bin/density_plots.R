#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rg",        type="character", default=NULL, help="Path to RGset RDS file"),
  make_option("--beta-norm", type="character", default=NULL, dest="beta_norm", help="Path to normalized beta RDS file")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path   <- opt$rg
beta_path <- opt$beta_norm
if (is.null(rg_path)   || !file.exists(rg_path))   stop("RGset RDS not found: ", rg_path)
if (is.null(beta_path) || !file.exists(beta_path)) stop("Beta norm RDS not found: ", beta_path)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

if (!requireNamespace("minfi", quietly=TRUE)) stop("Missing package: minfi")
suppressPackageStartupMessages(library(minfi))

tryCatch({
  rg <- readRDS(rg_path)
  beta_norm <- readRDS(beta_path)

  say("== Raw beta density (minfi) ==")
  beta_raw <- minfi::getBeta(minfi::preprocessRaw(rg))
  save_png("06_density_beta_raw.png",
           minfi::densityPlot(beta_raw, main="Raw beta distributions (after sample filtering)", xlab="Beta value"),
           w=1400, h=900)

  say("== Normalized beta density (SeSAMe) ==")
  save_png("07_density_beta_sesame.png",
           minfi::densityPlot(beta_norm, main="SeSAMe betas + pOOBAH QC", xlab="Beta value"),
           w=1400, h=900)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

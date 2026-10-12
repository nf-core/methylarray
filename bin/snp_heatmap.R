#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rg", type="character", default=NULL, help="Path to RGset RDS file")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path <- opt$rg
if (is.null(rg_path) || !file.exists(rg_path)) stop("RGset RDS not found: ", rg_path)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

if (!requireNamespace("minfi", quietly=TRUE)) stop("Missing package: minfi")
suppressPackageStartupMessages(library(minfi))

tryCatch({
  rg <- readRDS(rg_path)

  say("== SNP probes correlation heatmap ==")
  snp_beta <- try(minfi::getSnpBeta(rg), silent=TRUE)

  if (inherits(snp_beta, "try-error") || is.null(snp_beta) || !requireNamespace("pheatmap", quietly=TRUE)) {
    say("Skipping SNP heatmap (getSnpBeta and/or pheatmap not available).")
    quit(status=0)
  }

  cor_mat <- cor(snp_beta, use="pairwise.complete.obs")
  png("05_heatmap_snp_correlation.png", width=1200, height=1000, res=140)
  pheatmap::pheatmap(cor_mat, main="SNP probe correlation (sample identity check)")
  dev.off()
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

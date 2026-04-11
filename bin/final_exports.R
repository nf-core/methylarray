#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--beta-clean",   type="character", default=NULL,    dest="beta_clean",   help="Path to clean beta RDS file"),
  make_option("--export-csv",   type="character", default="false", dest="export_csv",   help="Export CSV files (true/false)"),
  make_option("--save-rds",     type="character", default="true",  dest="save_rds",     help="Save RDS files (true/false)"),
  make_option("--rds-compress", type="character", default="xz",    dest="rds_compress", help="RDS compression method")
)
opt <- parse_args(OptionParser(option_list=option_list))

beta_path  <- opt$beta_clean
export_csv <- tolower(opt$export_csv) %in% c("true","1","yes")
save_rds   <- tolower(opt$save_rds)   %in% c("true","1","yes")
compress   <- opt$rds_compress

if (is.null(beta_path) || !file.exists(beta_path)) stop("Clean beta RDS not found: ", beta_path)

suppressPackageStartupMessages(library(readr))

tryCatch({
  beta_clean <- readRDS(beta_path)

  bVals <- beta_clean
  bt <- 1e-6
  b_clip <- pmin(pmax(bVals, bt), 1 - bt)
  mVals <- log2(b_clip / (1 - b_clip))

  readr::write_csv(data.frame(Sample=colnames(bVals)), "samples_final.csv")
  readr::write_csv(data.frame(Probe=rownames(bVals)), "probes_final.csv")

  saveRDS(bVals, "bVals_final.rds", compress=compress)
  saveRDS(mVals, "mVals_final.rds", compress=compress)

  if (export_csv) {
    readr::write_csv(as.data.frame(mVals), "mVals_final.csv")
    readr::write_csv(as.data.frame(bVals), "bVals_final.csv")
  }

  qc <- c(
    paste("Final beta dim:", paste(dim(bVals), collapse=" x ")),
    "==== QC SUMMARY ====",
    paste("Final samples:", ncol(bVals)),
    paste("Final probes :", nrow(bVals))
  )
  writeLines(qc, "qc_summary.txt")
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

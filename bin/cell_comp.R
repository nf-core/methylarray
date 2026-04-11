#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rg",  type="character", default=NULL,                    help="Path to RGset RDS file"),
  make_option("--out", type="character", default="cell_counts_blood.csv", help="Output CSV for cell counts")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path <- opt$rg
out_csv <- opt$out

if (is.null(rg_path) || !file.exists(rg_path)) stop("RGset RDS not found: ", rg_path)

suppressPackageStartupMessages(library(readr))

tryCatch({
  if (!requireNamespace("FlowSorted.Blood.EPIC", quietly=TRUE) ||
      !requireNamespace("minfi", quietly=TRUE)) {
    stop("Required packages FlowSorted.Blood.EPIC and/or minfi are not available.")
  }

  suppressPackageStartupMessages(library(minfi))
  suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))

  # Force local hub to avoid network calls
  if (requireNamespace("AnnotationHub",  quietly=TRUE)) AnnotationHub::setAnnotationHubOption("LOCAL", TRUE)
  if (requireNamespace("ExperimentHub", quietly=TRUE)) ExperimentHub::setExperimentHubOption("LOCAL", TRUE)

  rg <- readRDS(rg_path)
  mset_cc <- minfi::preprocessNoob(rg)
  beta_cc <- minfi::getBeta(mset_cc)

  pfx <- sub("(_.*)$", "", rownames(beta_cc))
  sum_beta <- rowsum(beta_cc, group=pfx, reorder=FALSE)
  n_per_pfx <- as.numeric(table(pfx)); names(n_per_pfx) <- names(table(pfx))
  n_vec <- n_per_pfx[rownames(sum_beta)]
  beta_collapsed <- sweep(sum_beta, 1, n_vec, FUN="/")

  idol <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
  comp <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable

  idol2 <- intersect(idol, rownames(beta_collapsed))
  if (length(idol2) < 200) {
    stop("Insufficient overlap with IDOL CpGs (n=", length(idol2), "). Cannot estimate cell composition.")
  }

  beta_idol <- beta_collapsed[idol2, , drop=FALSE]
  comp_idol <- comp[idol2, , drop=FALSE]

  props <- FlowSorted.Blood.EPIC::projectCellType_CP(
    beta_idol, comp_idol,
    contrastWBC=NULL, nonnegative=TRUE, lessThanOne=FALSE
  )

  props_df <- as.data.frame(props)
  props_df$Sample <- rownames(props_df)
  readr::write_csv(props_df, out_csv)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

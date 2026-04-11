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

  say("== qcReport + control strip plots ==")
  try(minfi::qcReport(rg, pdf="00_qcReport.pdf"), silent=TRUE)

  controls_to_plot <- c(
    "STAINING","EXTENSION","HYBRIDIZATION","TARGET REMOVAL",
    "SPECIFICITY I","SPECIFICITY II","NON-POLYMORPHIC","NEGATIVE"
  )

  for (cc in controls_to_plot) {
    fname <- paste0("00_controlStrip_", gsub("[^A-Za-z0-9]+", "_", cc), ".png")
    save_png(fname, try(minfi::controlStripPlot(rg, controls=cc), silent=TRUE))
  }

  save_png("03_control_bisulfite_I.png",
           try(minfi::controlStripPlot(rg, controls="BISULFITE CONVERSION I"), silent=TRUE),
           w=1400, h=900, res=140)

  save_png("04_control_bisulfite_II.png",
           try(minfi::controlStripPlot(rg, controls="BISULFITE CONVERSION II"), silent=TRUE),
           w=1400, h=900, res=140)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

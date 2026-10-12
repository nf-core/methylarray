#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--sample-sheet", type="character", default=NULL,   dest="sample_sheet", help="Path to sample sheet TSV"),
  make_option("--idat-dir",     type="character", default=NULL,   dest="idat_dir",     help="Path to IDAT directory"),
  make_option("--annotation",   type="character", default="20a1.hg38", help="EPICv2 annotation string")
)
opt <- parse_args(OptionParser(option_list=option_list))

sample_sheet <- opt$sample_sheet
idat_dir     <- opt$idat_dir
annotation_v <- opt$annotation

if (is.null(sample_sheet) || !file.exists(sample_sheet)) stop("Sample sheet not found: ", sample_sheet)
if (is.null(idat_dir) || !dir.exists(idat_dir)) stop("IDAT dir not found: ", idat_dir)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

req("minfi")
req("readr")

force_epicv2_annotation <- function(rg, annotation_string) {
  try({
    annotation(rg) <- c(
      array = "IlluminaHumanMethylationEPICv2",
      annotation = annotation_string
    )
  }, silent = TRUE)
  rg
}

tryCatch({
  say("== Import IDATs + sample sheet ==")
  say("sample_sheet:", sample_sheet)
  say("idat_dir    :", idat_dir)

  ss <- readr::read_tsv(
    sample_sheet,
    col_types = readr::cols(
      Sample_Name      = readr::col_character(),
      Sentrix_ID       = readr::col_character(),
      Sentrix_Position = readr::col_character()
    ),
    trim_ws = TRUE
  )

  targets <- data.frame(
    Sample_Name      = ss$Sample_Name,
    Sentrix_ID       = trimws(ss$Sentrix_ID),
    Sentrix_Position = trimws(ss$Sentrix_Position),
    stringsAsFactors = FALSE
  )

  # Propagate optional sample_id column (explicit link to metadata)
  if ("sample_id" %in% names(ss)) {
    targets$Sample_id <- as.character(ss$sample_id)
    say("sample_id column found in samplesheet: will use it to join with metadata.")
  }

  targets$Basename <- file.path(
    idat_dir,
    targets$Sentrix_ID,
    paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)
  )

  missing_red <- !file.exists(paste0(targets$Basename, "_Red.idat"))
  missing_grn <- !file.exists(paste0(targets$Basename, "_Grn.idat"))
  if (any(missing_red | missing_grn)) {
    cat("Missing IDATs for:\n")
    print(targets[missing_red | missing_grn, c("Sample_Name", "Basename")])
    stop("Fix --idat_dir or the IDAT folder structure/names.")
  }

  dup_idx <- duplicated(targets$Basename) | duplicated(targets$Basename, fromLast = TRUE)
  if (any(dup_idx)) {
    cat("Duplicate basenames detected:\n")
    print(targets[dup_idx, c("Sample_Name", "Sentrix_ID", "Sentrix_Position", "Basename")])
    stop("Two rows map to the same Sentrix_ID + position.")
  }

  rg0 <- minfi::read.metharray.exp(targets = targets, extended = TRUE, verbose = TRUE)
  rg0 <- force_epicv2_annotation(rg0, annotation_v)
  minfi::sampleNames(rg0) <- targets$Sample_Name

  saveRDS(rg0, file = "rg0.rds")
  write.table(targets, file = "targets.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  say("Wrote: rg0.rds, targets.tsv")
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

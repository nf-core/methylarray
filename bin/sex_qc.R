#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rg",      type="character", default=NULL, help="Path to RGset RDS file"),
  make_option("--targets", type="character", default=NULL, help="Path to targets TSV"),
  make_option("--meta",    type="character", default="",   help="Path to metadata CSV (optional)")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path      <- opt$rg
targets_path <- opt$targets
meta_file    <- opt$meta

if (is.null(rg_path)      || !file.exists(rg_path))      stop("RGset RDS not found: ", rg_path)
if (is.null(targets_path) || !file.exists(targets_path)) stop("Targets TSV not found: ", targets_path)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

req("minfi"); req("readr"); req("dplyr")

tryCatch({
  rg <- readRDS(rg_path)
  targets <- readr::read_tsv(targets_path, show_col_types=FALSE)

  say("== Predicted sex (minfi::getSex) + optional metadata concordance ==")

  mset <- minfi::preprocessNoob(rg)
  if (!inherits(mset, c("GenomicMethylSet","GenomicRatioSet"))) mset <- minfi::mapToGenome(mset)

  pred_sex <- try(minfi::getSex(mset), silent=TRUE)
  targets$Sample_id <- targets$Sample_Name

  meta_tmp <- NULL
  if (nzchar(meta_file) && file.exists(meta_file)) {
    meta_tmp <- readr::read_csv(meta_file, show_col_types=FALSE) |> dplyr::rename_with(tolower)
  } else {
    say("Metadata not found or not provided -> skipping concordance.")
  }

  if (!inherits(pred_sex, "try-error")) {
    sex_df <- data.frame(
      Sample  = rownames(pred_sex),
      PredSex = pred_sex$predictedSex,
      xMed    = pred_sex$xMed,
      yMed    = pred_sex$yMed,
      stringsAsFactors = FALSE
    )
    sex_df$sample_id <- sex_df$Sample

    if (!is.null(meta_tmp) && all(c("sample_id","sex") %in% names(meta_tmp))) {
      meta_tmp$sample_id <- as.character(meta_tmp$sample_id)
      sex_df <- dplyr::left_join(sex_df, meta_tmp |> dplyr::select(sample_id, sex), by="sample_id")

      sex_df$MetaSex <- dplyr::case_when(
        sex_df$sex %in% c(0,"0","M","m","Male","male","MALE") ~ "M",
        sex_df$sex %in% c(1,"1","F","f","Female","female","FEMALE") ~ "F",
        TRUE ~ NA_character_
      )

      sex_df$Concordant <- ifelse(is.na(sex_df$MetaSex), NA, substr(sex_df$PredSex,1,1) == sex_df$MetaSex)
      discordant <- sex_df |> dplyr::filter(!is.na(MetaSex) & !Concordant)
      if (nrow(discordant) > 0) readr::write_csv(discordant, "sex_discordant_samples.csv")
    }

    readr::write_csv(sex_df, "qc_predicted_sex.csv")

    mset2 <- try(minfi::addSex(mset), silent=TRUE)
    if (!inherits(mset2, "try-error")) {
      save_png("08_plotSex_predicted.png", minfi::plotSex(mset2, id=minfi::sampleNames(mset2)))
    }
  } else {
    stop("getSex failed")
  }
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

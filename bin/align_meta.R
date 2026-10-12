#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--bvals",       type="character", default=NULL,                   help="Path to bVals RDS file"),
  make_option("--mvals",       type="character", default=NULL,                   help="Path to mVals RDS file"),
  make_option("--targets",     type="character", default=NULL,                   help="Path to targets TSV"),
  make_option("--meta",        type="character", default=NULL,                   help="Path to metadata CSV"),
  make_option("--cell-counts", type="character", default="",                     dest="cell_counts",  help="Path to cell counts CSV (optional)"),
  make_option("--pred-sex",    type="character", default="",                     dest="pred_sex",     help="Path to predicted sex CSV (optional)"),
  make_option("--keep-groups", type="character", default=NULL,                  dest="keep_groups",  help="Comma-separated list of groups to keep (default: keep all groups present in metadata)"),
  make_option("--export-csv",    type="character", default="false",     dest="export_csv",    help="Export aligned matrices as CSV (true/false)"),
  make_option("--include-array", type="character", default="false",     dest="include_array", help="Include Sentrix_Position (Array) as batch covariate (true/false)"),
  make_option("--cofactors",     type="character", default="sex,age,bmi", dest="cofactors",   help="Comma-separated cofactor names (case-insensitive). Standard: sex, age, bmi. Others must match metadata columns.")
)
opt <- parse_args(OptionParser(option_list=option_list))

bvals_path   <- opt$bvals
mvals_path   <- opt$mvals
targets_path <- opt$targets
meta_path    <- opt$meta
cell_path    <- opt$cell_counts
pred_path    <- opt$pred_sex
keep_groups  <- if (!is.null(opt$keep_groups) && nzchar(opt$keep_groups)) trimws(strsplit(opt$keep_groups, ",")[[1]]) else NULL
export_csv    <- tolower(opt$export_csv)    %in% c("true","1","yes")
include_array <- tolower(opt$include_array) %in% c("true","1","yes")

# Standard cofactor name → canonical column name produced by this script
COF_MAP       <- c(sex="Sex", age="Age", bmi="BMI")
cofactors_raw <- trimws(tolower(strsplit(opt$cofactors, ",")[[1]]))
cofactors_raw <- cofactors_raw[nzchar(cofactors_raw)]
extra_cofs    <- setdiff(cofactors_raw, names(COF_MAP))

if (is.null(bvals_path)   || !file.exists(bvals_path))   stop("bVals RDS not found: ", bvals_path)
if (is.null(mvals_path)   || !file.exists(mvals_path))   stop("mVals RDS not found: ", mvals_path)
if (is.null(targets_path) || !file.exists(targets_path)) stop("Targets TSV not found: ", targets_path)
if (is.null(meta_path)    || !file.exists(meta_path))    stop("Metadata CSV not found: ", meta_path)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

strip_leading_X_if_numeric <- function(x) sub("^X(?=\\d+$)", "", x, perl=TRUE)

cell_names <- c("CD8T","CD4T","NK","Bcell","Mono","Neu","Gran")

tryCatch({
  say("== Align matrices + metadata ==")

  bVals <- readRDS(bvals_path)
  mVals <- readRDS(mvals_path)
  targets <- readr::read_tsv(targets_path, show_col_types=FALSE)

  meta <- readr::read_csv(meta_path, show_col_types=FALSE) %>% dplyr::rename_with(tolower)

  # Always required: sample_id, group, sex
  hardcols     <- c("sample_id", "group", "sex")
  missing_hard <- setdiff(hardcols, names(meta))
  if (length(missing_hard) > 0) {
    stop("Missing required metadata columns: ", paste(missing_hard, collapse=", "))
  }

  # Cofactor columns: only required if explicitly requested via --cofactors
  cof_needed   <- cofactors_raw[cofactors_raw != "sex"]   # sex already checked above
  missing_cofs <- setdiff(cof_needed, names(meta))
  if (length(missing_cofs) > 0) {
    stop("Cofactor columns requested via --cofactors but missing from metadata: ",
         paste(missing_cofs, collapse=", "),
         ". Available columns: ", paste(names(meta), collapse=", "))
  }

  meta$sample_id <- as.character(meta$sample_id)
  if (any(duplicated(meta$sample_id))) stop("Duplicated sample_id in metadata.")

  if ("age" %in% names(meta)) meta$age <- suppressWarnings(as.numeric(meta$age))
  if ("bmi" %in% names(meta)) meta$bmi <- suppressWarnings(as.numeric(meta$bmi))

  if (!"Sample_id" %in% names(targets)) targets$Sample_id <- targets$Sample_Name
  targets$Sample_id <- as.character(targets$Sample_id)
  id_map <- setNames(as.character(targets$Sample_id), targets$Sample_Name)

  colnames(bVals) <- strip_leading_X_if_numeric(colnames(bVals))
  colnames(mVals) <- strip_leading_X_if_numeric(colnames(mVals))

  common_cols <- intersect(colnames(bVals), colnames(mVals))
  if (length(common_cols) == 0) stop("bVals and mVals share no common sample columns.")
  bVals <- bVals[, common_cols, drop=FALSE]
  mVals <- mVals[, common_cols, drop=FALSE]

  mapped_ids <- unname(id_map[common_cols])
  bad <- is.na(mapped_ids)

  if (anyDuplicated(mapped_ids[!is.na(mapped_ids)]) > 0) stop("Duplicated mapped sample_id after mapping.")

  if (any(bad)) {
    bad_names <- common_cols[bad]
    readr::write_csv(data.frame(Sample_Name=bad_names), "bVals_unmapped_colnames.csv")
    keep_cols <- common_cols[!bad]
    bVals <- bVals[, keep_cols, drop=FALSE]
    mVals <- mVals[, keep_cols, drop=FALSE]
    mapped_ids <- mapped_ids[!bad]
  }

  colnames(bVals) <- as.character(mapped_ids)
  colnames(mVals) <- as.character(mapped_ids)

  # Cell fractions: use only if file exists AND non-empty, otherwise ensure they are absent
  has_cell_file <- nzchar(cell_path) && file.exists(cell_path)

  if (has_cell_file) {
    props <- readr::read_csv(cell_path, show_col_types=FALSE)
    if (nrow(props) > 0) {
      if (!"Sample" %in% names(props)) names(props)[1] <- "Sample"
      props$sample_id <- as.character(unname(id_map[props$Sample]))
      props <- props %>% dplyr::filter(!is.na(sample_id))
      meta <- meta %>% dplyr::left_join(props %>% dplyr::select(-Sample), by="sample_id")
    } else {
      has_cell_file <- FALSE
    }
  }

  if (!has_cell_file) {
    # Remove any cell columns that may already exist in the input metadata
    meta <- meta %>% dplyr::select(-dplyr::any_of(cell_names))
  }

  # Groups
  meta$group <- toupper(trimws(meta$group))
  available_groups <- sort(unique(meta$group))
  say("Groups available in metadata: ", paste(available_groups, collapse=", "))
  if (!is.null(keep_groups)) {
    keep_groups <- toupper(trimws(keep_groups))
    missing_in_meta <- setdiff(keep_groups, available_groups)
    if (length(missing_in_meta) > 0)
      warning("Groups in --keep-groups not found in metadata: ", paste(missing_in_meta, collapse=", "))
    meta <- meta %>% dplyr::filter(group %in% keep_groups)
    say("Groups kept: ", paste(keep_groups, collapse=", "),
        " | Dropped: ", paste(setdiff(available_groups, keep_groups), collapse=", "))
  } else {
    say("--keep-groups not set: keeping all ", length(available_groups), " groups")
    keep_groups <- available_groups
  }

  common_samples <- intersect(colnames(bVals), meta$sample_id)
  if (length(common_samples) == 0) stop("No common samples between matrices and metadata after filtering.")

  meta <- meta[meta$sample_id %in% common_samples, , drop=FALSE]
  meta <- meta[match(common_samples, meta$sample_id), , drop=FALSE]

  bVals <- bVals[, meta$sample_id, drop=FALSE]
  mVals <- mVals[, meta$sample_id, drop=FALSE]

  # Standard cofactors — computed only when requested
  if ("sex" %in% cofactors_raw) {
    meta$Sex <- dplyr::case_when(
      meta$sex %in% c(0,"0","M","MALE","male")     ~ "M",
      meta$sex %in% c(1,"1","F","FEMALE","female") ~ "F",
      TRUE                                          ~ NA_character_
    )
    meta$Sex <- factor(meta$Sex, levels=c("F","M"))

    if (nzchar(pred_path) && file.exists(pred_path)) {
      pred <- readr::read_csv(pred_path, show_col_types=FALSE) %>% dplyr::rename_with(tolower)
      if (all(c("sample","predsex") %in% names(pred))) {
        pred$sample_id <- as.character(unname(id_map[pred$sample]))
        pred$Sex_pred  <- dplyr::case_when(
          tolower(pred$predsex) %in% c("male","m")   ~ "M",
          tolower(pred$predsex) %in% c("female","f") ~ "F",
          TRUE                                        ~ NA_character_
        )
        pred_map <- setNames(pred$Sex_pred, pred$sample_id)
        miss     <- is.na(meta$Sex)
        if (any(miss)) {
          meta$Sex[miss] <- pred_map[meta$sample_id[miss]]
          meta$Sex       <- factor(as.character(meta$Sex), levels=c("F","M"))
        }
      }
    }
  }

  if ("age" %in% cofactors_raw) meta$Age <- suppressWarnings(as.numeric(meta$age))
  if ("bmi" %in% cofactors_raw) meta$BMI <- suppressWarnings(as.numeric(meta$bmi))

  # Extra cofactors: validate presence in metadata, then pass through as-is
  if (length(extra_cofs) > 0) {
    missing_extra <- setdiff(extra_cofs, names(meta))
    if (length(missing_extra) > 0) {
      stop("Cofactor(s) not found in metadata: ", paste(missing_extra, collapse=", "),
           "\nAvailable columns: ", paste(names(meta), collapse=", "))
    }
    say(sprintf("Extra cofactors included: %s", paste(extra_cofs, collapse=", ")))
  }

  # Convert cell columns only if we actually have them
  cell_cols <- intersect(names(meta), cell_names)
  if (has_cell_file && length(cell_cols) > 0) {
    meta[cell_cols] <- lapply(meta[cell_cols], function(v) suppressWarnings(as.numeric(v)))
  }

  meta$group <- factor(meta$group, levels=keep_groups)

  # Batch: Plate (always) + Array (optional) from targets
  plate_cols <- c("sample_id", "Plate")
  if (include_array) plate_cols <- c(plate_cols, "Array")

  plate_df <- targets %>%
    dplyr::transmute(
      sample_id = as.character(Sample_id),
      Plate     = as.factor(Sentrix_ID),
      Array     = as.factor(Sentrix_Position)
    ) %>%
    dplyr::select(dplyr::all_of(plate_cols))

  meta <- dplyr::left_join(meta, plate_df, by="sample_id")

  readr::write_csv(meta, "meta_aligned.csv")
  saveRDS(bVals, "bVals_aligned.rds", compress="xz")
  saveRDS(mVals, "mVals_aligned.rds", compress="xz")

  if (export_csv) {
    readr::write_csv(tibble::as_tibble(bVals, rownames="Probe"), "bVals_aligned.csv")
    readr::write_csv(tibble::as_tibble(mVals, rownames="Probe"), "mVals_aligned.csv")
  }

  say("Aligned: meta_aligned.csv + bVals_aligned.rds + mVals_aligned.rds")
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--bvals",         type="character", default=NULL,          help="Path to bVals RDS file"),
  make_option("--meta",          type="character", default=NULL,          help="Path to metadata CSV"),
  make_option("--cofactors",     type="character", default="sex,age,bmi", help="Comma-separated cofactor names (same as passed to align_meta/dmp_limma)"),
  make_option("--n-probes-pca",  type="integer",   default=50000L,        dest="n_probes_pca", help="Number of most-variable probes for PCA"),
  make_option("--n-pcs",         type="integer",   default=10L,           dest="n_pcs",        help="Number of principal components to test"),
  make_option("--run-champ",     type="character", default="false",       dest="run_champ",    help="Run ChAMP SVD analysis (true/false)")
)
opt <- parse_args(OptionParser(option_list=option_list))

bvals_path  <- opt$bvals
meta_path   <- opt$meta
run_champ   <- tolower(opt$run_champ) %in% c("true","1","yes")
n_probes_pca <- as.integer(opt$n_probes_pca)
n_pcs        <- as.integer(opt$n_pcs)

# Standard cofactor name → canonical column name produced by align_meta.R
COF_MAP       <- c(sex="Sex", age="Age", bmi="BMI")
cofactors_raw <- trimws(tolower(strsplit(opt$cofactors, ",")[[1]]))
cofactors_raw <- cofactors_raw[nzchar(cofactors_raw)]
cof_cols      <- unname(ifelse(cofactors_raw %in% names(COF_MAP),
                               COF_MAP[cofactors_raw],
                               cofactors_raw))

if (is.null(bvals_path) || !file.exists(bvals_path)) stop("bVals RDS not found: ", bvals_path)
if (is.null(meta_path)  || !file.exists(meta_path))  stop("Metadata CSV not found: ", meta_path)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(matrixStats))

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

tryCatch({
  bVals <- readRDS(bvals_path)
  meta <- readr::read_csv(meta_path, show_col_types=FALSE)
  meta$sample_id <- as.character(meta$sample_id)

  if (!identical(colnames(bVals), as.character(meta$sample_id))) {
    stop("Column names of bVals do not match meta$sample_id.")
  }

  say("== Covariate checks + PCA association ==")

  cell_cols  <- intersect(names(meta), c("CD8T","CD4T","NK","Bcell","Mono","Neu","Gran"))
  std_num    <- intersect(cof_cols, c("Age","BMI"))
  extra_cofs <- setdiff(cof_cols, c("Sex","Age","BMI"))
  is_num_cof <- if (length(extra_cofs) > 0) vapply(extra_cofs, function(v) is.numeric(meta[[v]]), logical(1)) else logical(0)
  cov_num    <- c(std_num, cell_cols, extra_cofs[is_num_cof])
  cov_fac    <- c("group", "Plate",
                  if ("Sex" %in% cof_cols) "Sex" else character(0),
                  extra_cofs[!is_num_cof])

  na_tab <- sapply(c(cov_fac,cov_num), function(v) sum(is.na(meta[[v]])))
  txt <- capture.output({
    cat("-- NA counts --\n"); print(na_tab)
    cat("-- Group sizes --\n"); print(table(meta$group, useNA="ifany"))
    cat("-- Sex x group --\n"); print(table(meta$Sex, meta$group, useNA="ifany"))

    tmp <- meta %>% dplyr::filter(!is.na(Plate), !is.na(group))
    if (nrow(tmp) > 0) {
      cat("-- Plate x group (head) --\n")
      print(utils::head(table(tmp$Plate, tmp$group), 20))
      cat("[Plate ~ group] Chi-square test\n")
      print(chisq.test(table(tmp$Plate, tmp$group)))
    }

    if ("age" %in% names(meta)) {
      tmp <- meta %>% dplyr::filter(!is.na(age), !is.na(group))
      cat("[age ~ group] ANOVA\n"); print(summary(aov(age ~ group, data=tmp)))
    }
    if ("bmi" %in% names(meta)) {
      tmp <- meta %>% dplyr::filter(!is.na(bmi), !is.na(group))
      cat("[bmi ~ group] ANOVA\n"); print(summary(aov(bmi ~ group, data=tmp)))
    }
    if (length(cell_cols) > 0) {
      for (cc in cell_cols) {
        tmp <- meta %>% dplyr::filter(!is.na(.data[[cc]]), !is.na(group))
        cat("[", cc, " ~ group] ANOVA\n", sep="")
        print(summary(aov(tmp[[cc]] ~ tmp$group)))
      }
    }
  })
  writeLines(txt, "covariate_checks.txt")

  # PCA on top variable CpGs
  set.seed(1)
  n_select <- min(n_probes_pca, nrow(bVals))
  vars <- matrixStats::rowVars(bVals, na.rm=TRUE)
  sel <- order(vars, decreasing=TRUE)[seq_len(n_select)]

  beta_sub <- t(bVals[sel, , drop=FALSE])  # samples x probes
  cm <- colMeans(beta_sub, na.rm=TRUE)
  na_idx <- which(is.na(beta_sub), arr.ind=TRUE)
  if (nrow(na_idx) > 0) beta_sub[na_idx] <- cm[na_idx[,2]]

  pca <- prcomp(beta_sub, center=TRUE, scale.=FALSE)
  n_pcs_actual <- min(n_pcs, ncol(pca$x))
  pcs <- as.data.frame(pca$x[, seq_len(n_pcs_actual), drop=FALSE])
  pcs$sample_id <- rownames(pcs)
  meta_pcs <- dplyr::left_join(meta, pcs, by="sample_id")

  test_pc_assoc <- function(pc, varname, df) {
    x <- df[[varname]]; y <- df[[pc]]
    ok <- complete.cases(x,y)
    if (sum(ok) < 6) return(NA_real_)
    if (is.numeric(x)) {
      p <- try(summary(lm(y ~ x))$coefficients[2,4], silent=TRUE)
    } else {
      p <- try(anova(lm(y ~ as.factor(x)))$`Pr(>F)`[1], silent=TRUE)
    }
    if (inherits(p,"try-error")) return(NA_real_)
    as.numeric(p)
  }

  batch_cols   <- intersect(c("Plate","Array"), names(meta_pcs))
  extra_present <- intersect(extra_cofs, names(meta_pcs))
  vars_to_test <- unique(c("group", if ("Sex" %in% cof_cols) "Sex",
                            intersect(c("Age","BMI"), cof_cols),
                            cell_cols, batch_cols, extra_present))
  vars_to_test <- vars_to_test[vars_to_test %in% names(meta_pcs)]
  pc_names <- paste0("PC", seq_len(n_pcs_actual))

  pvals <- matrix(NA_real_, nrow=length(vars_to_test), ncol=length(pc_names),
                  dimnames=list(vars_to_test, pc_names))
  for (v in vars_to_test) for (pcn in pc_names) pvals[v, pcn] <- test_pc_assoc(pcn, v, meta_pcs)

  mlog <- -log10(pvals)
  write.csv(round(mlog,4), "covariates_PC_association_minuslog10p.csv")

  # Optional ChAMP SVD
  if (run_champ && requireNamespace("ChAMP", quietly=TRUE)) {
    out_dir <- "ChAMP_SVD"
    dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

    pd <- meta %>%
      dplyr::transmute(
        Sample_Name  = as.character(sample_id),
        Sample_Group = factor(group),
        Sex   = factor(Sex),
        Age   = suppressWarnings(as.numeric(age)),
        BMI   = suppressWarnings(as.numeric(bmi)),
        Plate = factor(Plate),
        dplyr::across(dplyr::all_of(intersect("Array", names(meta))), ~ factor(.x)),
        dplyr::across(dplyr::all_of(cell_cols), ~ suppressWarnings(as.numeric(.x)))
      )
    pd <- as.data.frame(pd); rownames(pd) <- pd$Sample_Name

    bVals2 <- bVals
    na_thr <- 0.01
    na_rate <- rowMeans(is.na(bVals2))
    bVals2 <- bVals2[na_rate <= na_thr, , drop=FALSE]
    if (anyNA(bVals2)) {
      rm <- rowMeans(bVals2, na.rm=TRUE)
      idx <- which(is.na(bVals2), arr.ind=TRUE)
      bVals2[idx] <- rm[idx[,1]]
    }

    ChAMP::champ.SVD(beta=as.data.frame(bVals2), pd=pd, resultsDir=out_dir)
  }
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

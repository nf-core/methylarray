#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--mvals",     type="character", default=NULL,          help="Path to mVals RDS file"),
  make_option("--bvals",     type="character", default=NULL,          help="Path to bVals RDS file"),
  make_option("--meta",      type="character", default=NULL,          help="Path to metadata CSV"),
  make_option("--cofactors", type="character", default="sex,age,bmi", help="Comma-separated cofactor names (must exist as columns in meta_aligned.csv)")
)
opt <- parse_args(OptionParser(option_list=option_list))

mvals_path <- opt$mvals
bvals_path <- opt$bvals
meta_path  <- opt$meta

# Standard cofactor name → canonical column name in meta_aligned.csv
COF_MAP       <- c(sex="Sex", age="Age", bmi="BMI")
cofactors_raw <- trimws(tolower(strsplit(opt$cofactors, ",")[[1]]))
cofactors_raw <- cofactors_raw[nzchar(cofactors_raw)]
cof_cols      <- unname(ifelse(cofactors_raw %in% names(COF_MAP),
                               COF_MAP[cofactors_raw],
                               cofactors_raw))

if (is.null(mvals_path) || !file.exists(mvals_path)) stop("mVals RDS not found: ", mvals_path)
if (is.null(bvals_path) || !file.exists(bvals_path)) stop("bVals RDS not found: ", bvals_path)
if (is.null(meta_path)  || !file.exists(meta_path))  stop("Metadata CSV not found: ", meta_path)

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(tibble))

# Génère toutes les paires de contrastes depuis un vecteur de groupes.
# Retourne une liste avec :
#   $matrix             — matrice de contrastes pour limma::contrasts.fit()
#   $contrast_to_groups — liste nommée : nom_contraste → c(g_num, g_denom)
make_all_pairwise_contrasts <- function(groups, design) {
  groups_safe  <- make.names(groups)
  idx_pairs    <- combn(seq_along(groups), 2, simplify = FALSE)
  contr_exprs  <- vapply(idx_pairs, function(idx)
    paste0("group", groups_safe[idx[2]], " - group", groups_safe[idx[1]]), character(1))
  contr_names  <- vapply(idx_pairs, function(idx)
    paste0(groups[idx[2]], "_vs_", groups[idx[1]]), character(1))
  mat          <- limma::makeContrasts(contrasts = contr_exprs, levels = design)
  colnames(mat) <- contr_names
  list(
    matrix             = mat,
    contrast_to_groups = setNames(
      lapply(idx_pairs, function(idx) c(groups[idx[2]], groups[idx[1]])),
      contr_names
    )
  )
}

tryCatch({
  mVals <- readRDS(mvals_path)
  bVals <- readRDS(bvals_path)
  meta <- readr::read_csv(
    meta_path,
    show_col_types = FALSE,
    col_types = readr::cols(sample_id = readr::col_character())
  )

  if (!identical(colnames(mVals), as.character(meta$sample_id))) {
    stop("Column names of mVals do not match meta$sample_id.")
  }
  if (!identical(colnames(bVals), as.character(meta$sample_id))) {
    stop("Column names of bVals do not match meta$sample_id.")
  }

  lambda_gc <- function(p){
    p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
    if (length(p) < 1000) return(NA_real_)
    stats::median(stats::qchisq(1 - p, df=1), na.rm=TRUE) / stats::qchisq(0.5, df=1)
  }

  run_dmp_limma <- function(method=c("plate_in_model","combat","norm_only")) {
    method <- match.arg(method)

    tryCatch({
      df <- meta
      df$group <- droplevels(factor(df$group))
      df$Plate <- droplevels(factor(df$Plate))

      # Coerce each cofactor to the appropriate type
      active_cofs <- cof_cols
      for (col in active_cofs) {
        if (!col %in% names(df)) stop("[DMP] Cofactor column not found in meta_aligned: ", col)
        if (col == "Sex") {
          df[[col]] <- droplevels(factor(df[[col]]))
        } else if (col %in% c("Age", "BMI")) {
          df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
        } else {
          num_try <- suppressWarnings(as.numeric(df[[col]]))
          if (sum(is.na(num_try) & !is.na(df[[col]])) == 0) {
            df[[col]] <- num_try
          } else {
            df[[col]] <- droplevels(factor(df[[col]]))
          }
        }
      }

      cell_all <- intersect(names(df), c("CD8T","CD4T","NK","Bcell","Mono","Neu","Gran"))
      cell_use <- cell_all
      if (length(cell_all) >= 2) {
        if ("Neu" %in% cell_all) cell_use <- setdiff(cell_all, "Neu") else cell_use <- cell_all[-length(cell_all)]
      }

      # Warn and drop single-level factors (would break full-rank)
      for (col in active_cofs) {
        if (is.factor(df[[col]]) && nlevels(df[[col]]) < 2) {
          warning(sprintf("[DMP] Cofactor '%s' has only 1 level — removing from design.", col))
          active_cofs <- setdiff(active_cofs, col)
        }
      }

      model_vars <- c("group", active_cofs, cell_use)
      if (method == "plate_in_model") model_vars <- c(model_vars, "Plate")

      keep_samples <- complete.cases(df[, model_vars, drop=FALSE])
      df2 <- droplevels(df[keep_samples, , drop=FALSE])

      mVals2 <- mVals[, df2$sample_id, drop=FALSE]
      bVals2 <- bVals[, df2$sample_id, drop=FALSE]

      na_thr <- 0.01
      bad <- !is.finite(mVals2) | is.na(mVals2)
      keep_cpg <- rowMeans(bad) <= na_thr
      mVals2 <- mVals2[keep_cpg, , drop=FALSE]
      bVals2 <- bVals2[rownames(mVals2), , drop=FALSE]

      if (any(!is.finite(mVals2) | is.na(mVals2))) {
        rm <- rowMeans(mVals2, na.rm=TRUE)
        idx <- which(!is.finite(mVals2) | is.na(mVals2), arr.ind=TRUE)
        mVals2[idx] <- rm[idx[,1]]
      }
      if (any(!is.finite(bVals2) | is.na(bVals2))) {
        rb <- rowMeans(bVals2, na.rm=TRUE)
        idx <- which(!is.finite(bVals2) | is.na(bVals2), arr.ind=TRUE)
        bVals2[idx] <- rb[idx[,1]]
      }

      rhs_base <- c("group", active_cofs)
      if (length(cell_use) > 0) rhs_base <- c(rhs_base, cell_use)

      rhs <- rhs_base
      if (method == "plate_in_model") rhs <- c(rhs, "Plate")

      out_dir <- paste0("DMP_limma__", method)
      dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

      # Cell-means coding (~ 0 + group): all group levels appear as named columns;
      # contrasts are built dynamically for all pairwise combinations of levels(df2$group).
      design <- model.matrix(as.formula(paste("~ 0 +", paste(rhs, collapse=" + "))), data=df2)
      colnames(design) <- make.names(colnames(design))
      if (!limma::is.fullrank(design)) stop("[DMP] Design not full-rank for method=", method)

      m_for_fit <- mVals2

      if (method == "combat") {
        if (!requireNamespace("sva", quietly=TRUE)) stop("[DMP] Package 'sva' required for ComBat method.")
        if (nlevels(df2$Plate) < 2) stop("[DMP] ComBat needs >=2 Plate levels.")

        mod_combat <- model.matrix(as.formula(paste("~ 0 +", paste(rhs_base, collapse=" + "))), data=df2)
        colnames(mod_combat) <- make.names(colnames(mod_combat))

        m_for_fit <- sva::ComBat(
          dat=as.matrix(mVals2),
          batch=df2$Plate,
          mod=mod_combat,
          par.prior=TRUE,
          prior.plots=FALSE
        )

        mod  <- mod_combat
        mod0 <- model.matrix(~ 1, data=df2)

        k_be   <- suppressWarnings(sva::num.sv(dat=m_for_fit, mod=mod, method="be"))
        k_leek <- suppressWarnings(sva::num.sv(dat=m_for_fit, mod=mod, method="leek"))
        k_est <- if (is.finite(k_be) && k_be > 0) k_be else if (is.finite(k_leek) && k_leek > 0) k_leek else 0L

        k_max <- max(0L, ncol(m_for_fit) - qr(mod)$rank - 1L)
        k_use <- min(as.integer(k_est), as.integer(k_max))

        design_try <- mod_combat
        if (k_use > 0) {
          svobj <- sva::sva(dat=m_for_fit, mod=mod, mod0=mod0, n.sv=k_use)
          sv <- svobj$sv; colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
          design_try <- cbind(mod_combat, sv)
        }

        while (k_use > 0 && !limma::is.fullrank(design_try)) {
          k_use <- k_use - 1L
          if (k_use == 0) { design_try <- mod_combat; break }
          svobj <- sva::sva(dat=m_for_fit, mod=mod, mod0=mod0, n.sv=k_use)
          sv <- svobj$sv; colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
          design_try <- cbind(mod_combat, sv)
        }

        design <- design_try
        if (!limma::is.fullrank(design)) stop("[DMP] Design not full-rank after SVA backoff.")

        write.csv(
          data.frame(n_samples=ncol(m_for_fit), mod_rank=qr(mod)$rank, k_be=k_be, k_leek=k_leek, k_est=k_est, k_max=k_max, k_used=k_use),
          file.path(out_dir, "SVA_k_selection.csv"),
          row.names=FALSE
        )
      }

      fit <- limma::lmFit(m_for_fit, design)

      pairwise           <- make_all_pairwise_contrasts(levels(df2$group), design)
      contr              <- pairwise$matrix
      contrast_to_groups <- pairwise$contrast_to_groups

      fit2 <- limma::contrasts.fit(fit, contr)
      fit2 <- limma::eBayes(fit2, trend=TRUE, robust=TRUE)

      get_delta_beta <- function(g1, g0) {
        rowMeans(bVals2[, df2$group == g1, drop=FALSE]) - rowMeans(bVals2[, df2$group == g0, drop=FALSE])
      }

      summ <- data.frame(contrast=colnames(contr), n_FDR_0.05=NA_integer_, lambda_gc=NA_real_, stringsAsFactors=FALSE)

      for (cn in colnames(contr)) {
        tryCatch({
          tt <- limma::topTable(fit2, coef=cn, number=Inf, adjust.method="BH", sort.by="P")
          g1 <- contrast_to_groups[[cn]][1]; g0 <- contrast_to_groups[[cn]][2]
          tt$deltaBeta <- get_delta_beta(g1, g0)[rownames(tt)]

          readr::write_csv(tibble::rownames_to_column(tt, "CpG"), file.path(out_dir, paste0("DMP_", cn, ".csv")))

          summ$n_FDR_0.05[summ$contrast==cn] <- sum(tt$adj.P.Val < 0.05, na.rm=TRUE)
          summ$lambda_gc[summ$contrast==cn]  <- lambda_gc(tt$P.Value)
        }, error = function(e) {
          message("ERROR in contrast ", cn, " (method=", method, "): ", conditionMessage(e))
          write.csv(data.frame(), file.path(out_dir, paste0("DMP_", cn, ".csv")), row.names=FALSE)
          summ$n_FDR_0.05[summ$contrast==cn] <<- 0L
        })
      }

      readr::write_csv(summ, file.path(out_dir, "DMP_summary_counts_and_lambda.csv"))
      readr::write_csv(df2, file.path(out_dir, "meta_used.csv"))
      summ
    }, error = function(e) {
      message("ERROR in run_dmp_limma(method=", method, "): ", conditionMessage(e))
      out_dir <- paste0("DMP_limma__", method)
      dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
      readr::write_csv(
        data.frame(contrast=character(), n_FDR_0.05=integer(), lambda_gc=double()),
        file.path(out_dir, "DMP_summary_counts_and_lambda.csv")
      )
    })
  }

  dir.create("DMP_limma__plate_in_model", showWarnings=FALSE, recursive=TRUE)
  dir.create("DMP_limma__combat", showWarnings=FALSE, recursive=TRUE)
  dir.create("DMP_limma__norm_only", showWarnings=FALSE, recursive=TRUE)

  run_dmp_limma("plate_in_model")
  run_dmp_limma("combat")
  run_dmp_limma("norm_only")
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

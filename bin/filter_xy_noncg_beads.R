#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

option_list <- list(
  make_option("--rg",             type="character", default=NULL,                    help="Path to RGChannelSet RDS (required when --filter-beads true)"),
  make_option("--beta-norm",      type="character", default=NULL,  dest="beta_norm", help="Path to normalized beta RDS file"),
  make_option("--annotation",     type="character", default="20a1.hg38",             help="EPICv2 annotation string"),
  make_option("--out-beta",       type="character", default="beta3.rds",             dest="out_beta",    help="Output filtered beta RDS"),
  make_option("--out-summary",    type="character", default="filtering_summary.tsv", dest="out_summary", help="Output filtering summary TSV"),
  make_option("--filter-xy",      type="character", default="true",  dest="filter_xy",      help="Filter XY probes (true/false)"),
  make_option("--filter-noncg",   type="character", default="true",  dest="filter_noncg",   help="Filter non-CpG probes (true/false)"),
  make_option("--filter-beads",   type="character", default="true",  dest="filter_beads",   help="Filter low-bead probes (true/false)"),
  make_option("--bead-min",       type="integer",   default=3L,      dest="bead_min",       help="Minimum beads per probe per sample"),
  make_option("--bead-fail-frac", type="double",    default=0.05,    dest="bead_fail_frac", help="Max fraction of samples below bead_min to retain probe")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg_path        <- opt$rg
beta_path      <- opt$beta_norm
anno_str       <- opt$annotation
out_beta       <- opt$out_beta
out_sum        <- opt$out_summary
filter_xy      <- tolower(opt$filter_xy)    %in% c("true","1","yes","y","t")
filter_noncg   <- tolower(opt$filter_noncg) %in% c("true","1","yes","y","t")
filter_beads   <- tolower(opt$filter_beads) %in% c("true","1","yes","y","t")
bead_min       <- opt$bead_min
bead_fail_frac <- opt$bead_fail_frac

if (filter_beads && (is.null(rg_path) || !file.exists(rg_path))) {
  stop("--rg must point to a valid RGChannelSet RDS when --filter-beads true")
}
if (is.null(beta_path) || !file.exists(beta_path)) stop("Beta norm RDS not found: ", beta_path)

req("minfi")
anno_pkg <- paste0("IlluminaHumanMethylationEPICv2anno.", anno_str)
req(anno_pkg)

tryCatch({
  beta_norm <- readRDS(beta_path)

  say("== Optional filters: beads=", filter_beads,
      " XY=", filter_xy,
      " nonCG=", filter_noncg, " ==")

  # ── 0) Bead count filter ────────────────────────────────────────────────────
  # Applied on the full beta matrix before probe-annotation intersection.
  # Uses minfi::getNBeads() on rg_filtered.rds (post-sample-QC RGChannelSet,
  # matching the samples in beta_norm). getNBeads() returns bead counts per
  # array address; addresses are mapped to probe names via getProbeInfo().
  # For Type I probes (two addresses per probe), the minimum bead count across
  # both the M and U addresses is used as the per-probe-per-sample count.
  # A probe is removed if > bead_fail_frac of samples fall below bead_min.
  beta_pre      <- beta_norm
  dropped_beads <- 0L
  if (isTRUE(filter_beads)) {
    rg <- readRDS(rg_path)

    # Align rg samples to beta_norm samples
    rg_samps  <- colnames(rg)
    bet_samps <- colnames(beta_pre)
    common_s  <- intersect(rg_samps, bet_samps)
    if (length(common_s) == 0) {
      say("WARNING: no common samples between RGset and beta_norm -> skipping bead filter")
    } else {
      rg <- rg[, common_s]

      nb <- try(minfi::getNBeads(rg), silent=TRUE)
      if (!inherits(nb, "try-error") && !is.null(nb) && nrow(nb) > 0) {
        nb_addr <- rownames(nb)

        # Type II probes: one address (AddressA) per probe
        ti2    <- minfi::getProbeInfo(rg, type="II")
        t2_map <- setNames(as.character(ti2$Name), as.character(ti2$AddressA))
        t2_addr <- intersect(nb_addr, names(t2_map))
        nb_t2   <- nb[t2_addr, , drop=FALSE]
        rownames(nb_t2) <- t2_map[t2_addr]

        # Type I probes: two addresses (AddressA = M channel, AddressB = U channel)
        # Use min bead count across both channels per probe
        ti1    <- minfi::getProbeInfo(rg, type="I")
        t1_ma  <- setNames(as.character(ti1$Name), as.character(ti1$AddressA))
        t1_ub  <- setNames(as.character(ti1$Name), as.character(ti1$AddressB))
        t1a_addr <- intersect(nb_addr, names(t1_ma))
        t1b_addr <- intersect(nb_addr, names(t1_ub))
        if (length(t1a_addr) > 0 && length(t1b_addr) > 0) {
          nb_t1a <- nb[t1a_addr, , drop=FALSE]; rownames(nb_t1a) <- t1_ma[t1a_addr]
          nb_t1b <- nb[t1b_addr, , drop=FALSE]; rownames(nb_t1b) <- t1_ub[t1b_addr]
          common_t1 <- intersect(rownames(nb_t1a), rownames(nb_t1b))
          nb_t1 <- pmin(nb_t1a[common_t1, , drop=FALSE], nb_t1b[common_t1, , drop=FALSE])
        } else {
          nb_t1 <- nb[character(0), , drop=FALSE]
        }

        nb_by_probe <- rbind(nb_t2, nb_t1)
        common_nb   <- intersect(rownames(nb_by_probe), rownames(beta_pre))

        if (length(common_nb) > 0) {
          nb_sub <- nb_by_probe[common_nb, common_s, drop=FALSE]
          nb_sub[is.na(nb_sub)] <- 0L
          low_beads <- rowMeans(nb_sub < bead_min) > bead_fail_frac
          drop_ids  <- common_nb[low_beads]
          dropped_beads <- length(drop_ids)
          beta_pre <- beta_pre[!(rownames(beta_pre) %in% drop_ids), , drop=FALSE]
          say("Bead filter: dropped", dropped_beads, "probes (bead_min=", bead_min,
              ", fail_frac>", bead_fail_frac, ")")
        } else {
          say("WARNING: no probe-address overlap found -> skipping bead filter")
        }
      } else {
        say("WARNING: getNBeads() unavailable -> skipping bead filter")
      }
    }
  }

  anno <- minfi::getAnnotation(get(anno_pkg))

  # ── 1) Restrict to probes in annotation ─────────────────────────────────────
  present <- intersect(rownames(anno), rownames(beta_pre))
  beta0   <- beta_pre[present, , drop=FALSE]

  # ── 2) XY chromosomes ───────────────────────────────────────────────────────
  beta1      <- beta0
  dropped_xy <- 0L
  if (isTRUE(filter_xy)) {
    xy_ids     <- present[anno[present, "chr"] %in% c("chrX","chrY")]
    dropped_xy <- length(intersect(rownames(beta1), xy_ids))
    beta1      <- beta1[!(rownames(beta1) %in% xy_ids), , drop=FALSE]
  }

  # ── 3) Non-CpG probes ───────────────────────────────────────────────────────
  beta2        <- beta1
  dropped_noncg <- 0L
  if (isTRUE(filter_noncg)) {
    cg_ids        <- rownames(beta2)[grepl("^cg", rownames(beta2))]
    dropped_noncg <- nrow(beta2) - length(cg_ids)
    beta2         <- beta2[cg_ids, , drop=FALSE]
  }

  beta3 <- beta2

  # ── Summary ─────────────────────────────────────────────────────────────────
  summ <- data.frame(
    stage = c(
      "input_beta_norm",
      if (filter_beads) "after_drop_beads" else "skip_drop_beads",
      "after_present",
      if (filter_xy)    "after_drop_xy"    else "skip_drop_xy",
      if (filter_noncg) "after_drop_non_cg" else "skip_drop_non_cg",
      "final_beta3"
    ),
    n_probes = c(
      nrow(beta_norm),
      nrow(beta_pre),
      nrow(beta0),
      nrow(beta1),
      nrow(beta2),
      nrow(beta3)
    ),
    n_samples     = rep(ncol(beta_norm), 6),
    dropped_beads = c(NA, dropped_beads, NA, NA, NA, NA),
    dropped_xy    = c(NA, NA, NA, dropped_xy, NA, NA),
    dropped_noncg = c(NA, NA, NA, NA, dropped_noncg, NA),
    stringsAsFactors = FALSE
  )

  write.table(summ, out_sum, sep="\t", quote=FALSE, row.names=FALSE)
  saveRDS(beta3, out_beta)

  say("Done: wrote ", out_beta, " and ", out_sum)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

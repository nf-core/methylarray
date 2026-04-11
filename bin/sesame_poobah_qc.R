#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--rg0",             type="character", default=NULL,          help="Path to RG0 RDS file"),
  make_option("--targets",         type="character", default=NULL,          help="Path to targets TSV"),
  make_option("--sesame-prep",     type="character", default="QCDPB",       dest="sesame_prep",    help="SeSAMe prep pipeline string"),
  make_option("--p-det",           type="double",    default=0.05,          dest="p_det",          help="pOOBAH p-value threshold"),
  make_option("--fail-sample-fr",  type="double",    default=0.10,          dest="fail_sample_fr", help="Max fraction of failed probes per sample"),
  make_option("--fail-probe-fr",   type="double",    default=0.10,          dest="fail_probe_fr",  help="Max fraction of failed samples per probe"),
  make_option("--out-beta",        type="character", default="beta_norm.rds",       dest="out_beta",    help="Output beta RDS"),
  make_option("--out-rg",          type="character", default="rg_filtered.rds",     dest="out_rg",      help="Output filtered RG RDS"),
  make_option("--out-targets",     type="character", default="targets_filtered.tsv", dest="out_targets", help="Output filtered targets TSV")
)
opt <- parse_args(OptionParser(option_list=option_list))

rg0_path     <- opt$rg0
targets_path <- opt$targets
sesame_prep  <- opt$sesame_prep
p_det        <- opt$p_det
fail_s_fr    <- opt$fail_sample_fr
fail_p_fr    <- opt$fail_probe_fr
out_beta     <- opt$out_beta
out_rg       <- opt$out_rg
out_targets  <- opt$out_targets

if (!file.exists(rg0_path))     stop("RG0 RDS not found: ", rg0_path)
if (!file.exists(targets_path)) stop("Targets TSV not found: ", targets_path)

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

req("readr"); req("minfi"); req("sesame"); req("sesameData")

# SeSAMe helper — applies the prep pipeline via openSesame().
# The Q step uses SeSAMe's pre-computed quality masks (mask_general etc.).
# Per-sample bead count filtering is handled downstream in FILTER_XY_NONCG_BEADS
# via minfi::getNBeads() on the filtered RGChannelSet.
sesame_betas <- function(idats_named, prep="QCDPB", pval_threshold=0.05){
  args <- list(
    x = idats_named,
    prep = prep,
    prep_args = list(
      P = list(pval.threshold = pval_threshold)
    ),
    func = sesame::getBetas
  )
  if (requireNamespace("BiocParallel", quietly=TRUE)) {
    args$BPPARAM <- BiocParallel::SerialParam()
  }
  do.call(sesame::openSesame, args)
}

tryCatch({
  say(paste0("== SeSAMe betas (", sesame_prep, ") + masking(P=", p_det, ") =="))
  try(sesameData::sesameDataCache(), silent=TRUE)

  rg0 <- readRDS(rg0_path)
  targets <- readr::read_tsv(targets_path, show_col_types=FALSE)

  # Build basenames from targets.tsv (must contain Basename column)
  if (!("Basename" %in% names(targets))) stop("targets.tsv must contain a 'Basename' column")

  idats <- targets$Basename
  names(idats) <- targets$Sample_Name

  beta_ses <- sesame_betas(idats, prep=sesame_prep, pval_threshold=p_det)

  if (is.null(colnames(beta_ses)) || any(colnames(beta_ses)=="")) colnames(beta_ses) <- names(idats)
  beta_ses <- beta_ses[, targets$Sample_Name, drop=FALSE]

  # Sample filtering
  fail_rate_sample <- colMeans(is.na(beta_ses))
  keep_samp <- fail_rate_sample <= fail_s_fr
  say("Samples: initial =", ncol(beta_ses), "| kept =", sum(keep_samp))
  if (any(!keep_samp)) {
    dropped <- data.frame(Sample=names(keep_samp)[!keep_samp], FailRate=fail_rate_sample[!keep_samp])
    readr::write_csv(dropped, "dropped_samples.csv")
    say("Dropped samples:", paste(dropped$Sample, collapse=", "))
  } else {
    readr::write_csv(data.frame(Sample=character(), FailRate=numeric()), "dropped_samples.csv")
  }

  beta1 <- beta_ses[, keep_samp, drop=FALSE]

  # Probe filtering
  fail_rate_probe <- rowMeans(is.na(beta1))
  keep_probe <- fail_rate_probe <= fail_p_fr
  beta_norm <- beta1[keep_probe, , drop=FALSE]

  say("Probes: initial =", nrow(beta_ses),
      "| after sample filter =", nrow(beta1),
      "| after probe filter =", nrow(beta_norm))

  # QC outputs
  save_png("02_hist_failrate_per_sample_masking.png", {
    hist(fail_rate_sample, breaks=30,
         main="Per-sample failed probe rate (SeSAMe masking)",
         xlab="Proportion of probes masked as NA in beta matrix")
    abline(v=fail_s_fr, lty=2)
  })

  readr::write_csv(
    data.frame(Sample=names(fail_rate_sample), FailRate=as.numeric(fail_rate_sample)),
    "masking_failrate_per_sample.csv"
  )
  readr::write_csv(
    data.frame(Probe=rownames(beta1), FailRate=as.numeric(fail_rate_probe)),
    "masking_failrate_per_probe.csv"
  )

  # Apply sample filtering to rg + targets
  rg <- rg0[, keep_samp]
  targets_f <- targets[keep_samp, , drop=FALSE]

  saveRDS(beta_norm, out_beta)
  saveRDS(rg, out_rg)
  write.table(targets_f, out_targets, sep="\t", quote=FALSE, row.names=FALSE)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

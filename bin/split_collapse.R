#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

.bin_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[[1]])))
source(file.path(.bin_dir, "utils.R"))

option_list <- list(
  make_option("--beta3",          type="character", default=NULL,                          help="Path to beta3 RDS file"),
  make_option("--out-dmr",        type="character", default="beta_for_DMR_uncollapsed.rds", dest="out_dmr",   help="Output DMR beta RDS (uncollapsed, no SNP/CH filtering)"),
  make_option("--out-clean",      type="character", default="beta_clean.rds",              dest="out_clean", help="Output clean DMP beta RDS (SNP/CH filtered + collapsed)"),
  make_option("--do-rmsnpandch",  type="character", default="true",  dest="do_rmsnpandch",  help="Run rmSNPandCH (true/false)"),
  make_option("--snp-dist-bp",    type="double",    default=2,       dest="snp_dist_bp",    help="SNP distance in bp for rmSNPandCH"),
  make_option("--snp-maf",        type="double",    default=0.01,    dest="snp_maf",        help="MAF cutoff for rmSNPandCH")
)
opt <- parse_args(OptionParser(option_list=option_list))

beta_path     <- opt$beta3
out_dmr       <- opt$out_dmr
out_clean     <- opt$out_clean
do_rmsnpandch <- tolower(opt$do_rmsnpandch) %in% c("true","1","yes")
snp_dist_bp   <- opt$snp_dist_bp
snp_maf       <- opt$snp_maf

if (is.null(beta_path) || !file.exists(beta_path)) stop("Beta3 RDS not found: ", beta_path)

req("DMRcate"); req("sesame")

tryCatch({
  beta_in <- readRDS(beta_path)

  say("== Split pipelines: DMP vs DMR + collapse EPICv2 replicates ==")

  # DMR branch: SNP filtering without cross-hybridisation removal
  beta_dmr <- beta_in
  if (do_rmsnpandch) {
    beta_dmr <- DMRcate::rmSNPandCH(beta_dmr, dist=snp_dist_bp, mafcut=snp_maf, rmcrosshyb=FALSE)
  }
  saveRDS(beta_dmr, out_dmr)
  say("DMR beta saved (uncollapsed):", nrow(beta_dmr))

  # DMP branch: SNP filtering with cross-hybridisation removal
  beta_dmp <- beta_in
  if (do_rmsnpandch) {
    beta_dmp <- DMRcate::rmSNPandCH(beta_dmp, dist=snp_dist_bp, mafcut=snp_maf, rmcrosshyb=TRUE)
  }
  say("DMP after rmSNPandCH:", nrow(beta_dmp))

  beta_clean <- sesame::betasCollapseToPfx(beta_dmp)
  say("DMP after collapse:", nrow(beta_clean))

  saveRDS(beta_clean, out_clean)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

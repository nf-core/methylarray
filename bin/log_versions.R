#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--out-versions", type="character", default="package_versions.csv", dest="out_versions", help="Output CSV for package versions"),
  make_option("--out-session",  type="character", default="sessionInfo.txt",       dest="out_session",  help="Output file for sessionInfo")
)
opt <- parse_args(OptionParser(option_list=option_list))

out_versions <- opt$out_versions
out_session  <- opt$out_session

tryCatch({
  pkg_versions <- function(pkgs){
    do.call(rbind, lapply(pkgs, function(p){
      v <- tryCatch(as.character(utils::packageVersion(p)), error=function(e) NA_character_)
      data.frame(package=p, version=v, stringsAsFactors=FALSE)
    }))
  }

  bioc_ver <- NA_character_
  if (requireNamespace("BiocManager", quietly=TRUE)) {
    bioc_ver <- tryCatch(as.character(BiocManager::version()), error=function(e) NA_character_)
  }

  versions <- rbind(
    data.frame(package="R", version=as.character(getRversion()), stringsAsFactors=FALSE),
    data.frame(package="Bioconductor", version=bioc_ver, stringsAsFactors=FALSE),
    pkg_versions(c(
      "minfi","sesame","sesameData","DMRcate","limma","sva",
      "IlluminaHumanMethylationEPICv2manifest",
      "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
      "readr","dplyr","ggplot2","matrixStats","tibble",
      "FlowSorted.Blood.EPIC","ChAMP","pheatmap","missmethyl"
    ))
  )

  write.csv(versions, out_versions, row.names=FALSE)
  writeLines(capture.output(sessionInfo()), out_session)
}, error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})

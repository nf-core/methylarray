#!/usr/bin/env Rscript
# Shared utility functions sourced by all pipeline R scripts.

say <- function(...) message(paste(..., collapse = " "))

req <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

save_png <- function(path, expr, w = 1200, h = 900, res = 150) {
    png(path, w, h, res = res)
    on.exit(dev.off(), add = TRUE)
    force(expr)
}

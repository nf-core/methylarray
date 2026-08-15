# nf-core/methylarray — Docker image
# Base: Bioconductor 3.19 (R 4.4)
# Packages: minfi, sesame, sesameData, limma, DMRcate, FlowSorted.Blood.EPIC, ChAMP, coMethDMR

FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL org.opencontainers.image.authors="nf-core/methylarray contributors"
LABEL org.opencontainers.image.source="https://github.com/nf-core/methylarray"
LABEL org.opencontainers.image.description="R/Bioconductor environment for EPICv2 methylation array analysis"

# System dependencies for R packages with compiled components
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libhdf5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages (pinned to Bioc 3.19 / R 4.4)
RUN R -e " \
    options(repos = BiocManager::repositories()); \
    pkgs <- c( \
        'minfi', \
        'sesame', \
        'sesameData', \
        'limma', \
        'DMRcate', \
        'FlowSorted.Blood.EPIC', \
        'IlluminaHumanMethylationEPICv2anno.20a1.hg38', \
        'IlluminaHumanMethylationEPICv2manifest', \
        'ChAMP', \
        'sva', \
        'missMethyl', \
        'pheatmap', \
        'ggplot2', \
        'dplyr', \
        'readr', \
        'optparse', \
        'yaml' \
    ); \
    BiocManager::install(pkgs, ask = FALSE, update = FALSE); \
    "

# coMethDMR is used by DMR_DMRCATE when --run-cometh is enabled (default: true).
# Installed as its own layer so it doesn't invalidate the cache of the package
# block above when new dependencies are added here in the future.
RUN R -e " \
    options(repos = BiocManager::repositories()); \
    BiocManager::install('coMethDMR', ask = FALSE, update = FALSE); \
    "

# Verify key packages load correctly
RUN R -e " \
    library(minfi); \
    library(sesame); \
    library(limma); \
    library(DMRcate); \
    cat('All packages loaded successfully\n'); \
    "

# Set standard nf-core R environment variables
ENV R_PROFILE_USER="/.Rprofile"
ENV R_ENVIRON_USER="/.Renviron"

# Nextflow's docker profile runs the container as the host UID/GID
# (docker.runOptions = '-u $(id -u):$(id -g)'), which has no matching user
# entry in the image. Without a writable $HOME, ExperimentHub/AnnotationHub/
# sesameData fail to create their cache directory ("Permission denied").
# A world-writable HOME keeps this working for any UID.
ENV HOME=/home/methylarray
RUN mkdir -p /home/methylarray && chmod 777 /home/methylarray

CMD ["R"]

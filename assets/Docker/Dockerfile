FROM r-base:4.3.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libicu-dev \
    build-essential \
    libfontconfig1-dev \
    libfreetype6-dev \
    curl \
    wget \
    git \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

RUN R -e "install.packages(c('pkgdown', 'systemfonts', 'textshaping', 'ragg'), repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install(c('minfi', 'IlluminaHumanMethylationEPICmanifest', 'limma', 'missMethyl', 'IlluminaHumanMethylation450kmanifest', 'minfiData', 'DMRcate', 'dplyr', 'readr'))"

# Install preprocessCore to patch preprocessCore issues
RUN R -e "BiocManager::install('preprocessCore', configure.args = c(preprocessCore = '--disable-threading'), force=TRUE, update=TRUE, type='source')"

RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19.masked', 'BSgenome.Hsapiens.UCSC.hg38.masked', 'wateRmelon'))"

RUN git clone https://github.com/pjhop/DNAmCrosshyb.git /DNAmCrosshyb
RUN R CMD INSTALL /DNAmCrosshyb

RUN R -e "BiocManager::install(c('IlluminaHumanMethylationEPICanno.ilm10b2.hg19', 'FlowSorted.Blood.EPIC'))"
RUN R -e "BiocManager::install(c('FlowSorted.Blood.450k'))"

RUN R -e "BiocManager::install(c('ChAMP'))"

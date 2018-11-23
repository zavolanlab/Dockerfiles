##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="4"
LABEL software="R_DGE_DTU"
LABEL software.version="3.5.1"
LABEL software.description="DGE & DTU using different R libraries."
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license=""
LABEL software.tags="Transcriptomics"
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.5.1

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y --allow-downgrades build-essential curl zlib1g-dev software-properties-common libcurl4=7.61.0-1 libcurl4-openssl-dev libxml2-dev libssl-dev apt-transport-https libgfortran3 libopenblas-base

# apt-get build-dep libcurl4-openssl-dev

##### INSTALL R packages #####
RUN Rscript -e 'install.packages(c("optparse", "devtools", "RCurl", "BiocManager"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'BiocManager::install(c("DRIMSeq", "DEXSeq", "DESeq", "edgeR", "stageR", "tximport", "GenomicFeatures"), ask=FALSE, update = FALSE, version="3.8");'

#  && apt-get autoremove -y \
#  && apt-get clean \
#  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

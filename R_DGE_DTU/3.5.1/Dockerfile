##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="3"
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
  && apt-get install -y curl zlib1g-dev software-properties-common libcurl4 libxml2-dev libssl-dev apt-transport-https libgfortran3 libopenblas-base \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
  && Rscript -e 'install.packages(c("optparse", "gplots", "RColorBrewer", "statmod"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("rtracklayer", "edgeR", "DRIMSeq", "DESeq", "stageR", "DEXSeq", "tximport","GenomicFeatures"), ask=FALSE);'

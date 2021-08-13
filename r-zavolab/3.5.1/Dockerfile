##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="3"
LABEL software="3.5.1"
LABEL software.version="3.5.1"
LABEL software.description="R is a system for statistical computation and graphics."
LABEL software.website="https://www.r-project.org/"
LABEL software.documentation="https://www.r-project.org/"
LABEL software.license=""
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.5.1

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y curl zlib1g-dev software-properties-common libcurl4-gnutls-dev libxml2-dev libssl-dev apt-transport-https \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
  && Rscript -e 'install.packages(c("devtools", "optparse", "gplots", "RColorBrewer", "statmod", "rjson", "dplyr"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("rtracklayer", "edgeR"), ask=FALSE);'

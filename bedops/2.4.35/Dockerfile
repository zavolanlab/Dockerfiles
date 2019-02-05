##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="bedops"
LABEL software.version="2.4.35"
LABEL software.description="BEDOPS: high-performance genomic feature operations"
LABEL software.website="https://github.com/bedops/bedops"
LABEL software.documentation="https://bedops.readthedocs.io/"
LABEL software.license="GNU General Public License v2.0"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.4.35

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget \
  && wget https://github.com/bedops/bedops/releases/download/v${SOFTWARE_VERSION}/bedops_linux_x86_64-v${SOFTWARE_VERSION}.tar.bz2 \
  && tar jxvf bedops_linux_x86_64-v${SOFTWARE_VERSION}.tar.bz2 \
  && rm bedops_linux_x86_64-v${SOFTWARE_VERSION}.tar.bz2 \
  && apt-get autoremove -y wget \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

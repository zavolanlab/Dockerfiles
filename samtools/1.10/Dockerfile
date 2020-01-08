##### BASE IMAGE #####
FROM ubuntu:18.04


##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="samtools"
LABEL software.version="1.10"
LABEL software.description="Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format"
LABEL software.website="https://github.com/samtools/samtools"
LABEL software.documentation="https://github.com/samtools/samtools"
LABEL software.license="MIT"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.10

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev ca-certificates zlib1g-dev curl libcurl4-openssl-dev unzip autoconf \
  && apt-get clean \
  && curl -k -L https://github.com/samtools/samtools/releases/download/${SOFTWARE_VERSION}/samtools-${SOFTWARE_VERSION}.tar.bz2 -o samtools-${SOFTWARE_VERSION}.tar.bz2 \
  && tar xvjf samtools-${SOFTWARE_VERSION}.tar.bz2 \
  && cd samtools-${SOFTWARE_VERSION} \
  && make \
  && make install \
  && cd ../ \
  && cp samtools-${SOFTWARE_VERSION}/samtools /usr/local/bin/ \
  && rm -r samtools-${SOFTWARE_VERSION}.tar.bz2 \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

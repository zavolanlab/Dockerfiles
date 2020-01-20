##### BASE IMAGE #####
FROM ubuntu:18.04


##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="bedtools"
LABEL software.version="2.28.0"
LABEL software.description="The bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks"
LABEL software.website="http://bedtools.readthedocs.io/en/latest/"
LABEL software.documentation="http://quinlanlab.org/tutorials/bedtools/bedtools.html"
LABEL software.license="GNU Public License (Version 2)"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.28.0

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev ca-certificates zlib1g-dev curl unzip autoconf \
  && apt-get clean \
  && cd $HOME \
  && curl -k -L https://github.com/arq5x/bedtools2/archive/v${SOFTWARE_VERSION}.tar.gz -o bedtools2-v${SOFTWARE_VERSION}.tar.gz \
  && tar xzvf bedtools2-v${SOFTWARE_VERSION}.tar.gz \
  && rm bedtools2-v${SOFTWARE_VERSION}.tar.gz \
  && cd bedtools2-${SOFTWARE_VERSION} \
  && make \
  && make install \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

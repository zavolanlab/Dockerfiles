##### BASE IMAGE #####
FROM ubuntu:18.04


##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="salmon"
LABEL software.version="0.11.0"
LABEL software.description="Salmon: Fast, accurate and bias-aware transcript quantification from RNA-seq data"
LABEL software.website="https://combine-lab.github.io/salmon/getting_started/"
LABEL software.documentation="http://salmon.readthedocs.io/en/latest/"
LABEL software.license="GNU General Public License v3.0"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 0.11.0

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev ca-certificates zlib1g-dev curl unzip autoconf \
  && apt-get clean \
  && curl -k -L https://github.com/COMBINE-lab/salmon/archive/v${SALMON_VERSION}.tar.gz -o salmon-v${SALMON_VERSION}.tar.gz \
  && tar xzf salmon-v${SOFTWARE_VERSION}.tar.gz \
  && cd salmon-${SOFTWARE_VERSION} \
  && mkdir build \
  && cd build \
  && cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install \
  && rm -r salmon-v${SOFTWARE_VERSION}.tar.gz \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

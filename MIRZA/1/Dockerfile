##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="3"
LABEL software="MIRZA"
LABEL software.version="1"
LABEL software.description="A biophysical model of miRNA-target interaction and infer its energy parameters from Ago-CLIP data."
LABEL software.website="http://mirz.unibas.ch"
LABEL software.documentation=""
LABEL software.license="SPDX:GPL-3.0-or-later"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="jaroslaw.surkont@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="SPDX:GPL-3.0"

COPY Makefile /tmp/mirza/
COPY src /tmp/mirza/src

RUN apt-get update && \
  apt-get install --yes \
    make \
    g++ && \
  cd /tmp/mirza/ && \
  make && \
  mv MIRZA /usr/bin/ && \
  rm -rf /tmp/mirza && \
  apt-get remove --purge --yes \
    make \
    g++ && \
  apt-get autoremove --purge --yes

##### BASE IMAGE #####
FROM python:2.7-slim

##### METADATA #####
LABEL base.image="python:2.7-slim"
LABEL version="5"
LABEL software="rmats"
LABEL software.version="4.0.2"
LABEL software.description="MATS is a computational tool to detect differential alternative splicing events from RNA-Seq data."
LABEL software.website="http://rnaseq-mats.sourceforge.net/"
LABEL software.documentation="http://rnaseq-mats.sourceforge.net/user_guide.htm"
LABEL software.license=""
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 4.0.2

RUN apt-get update -y \
  && apt-get install -y wget libblas-dev liblapack-dev libgsl2 libgsl-dev gfortran gcc make \
  && pip2 install numpy \
  && wget https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS.${SOFTWARE_VERSION}.tgz/download \
  && tar -xzf download \
  && chmod +x /rMATS.${SOFTWARE_VERSION}/rMATS-turbo-Linux-UCS4/rmats.py \
  && chmod +x /rMATS.${SOFTWARE_VERSION}/rMATS-turbo-Linux-UCS2/rmats.py \
  && rm download \
  && cd /usr/lib/x86_64-linux-gnu \
  && ln -s libgsl.so libgsl.so.0 \
  && cd $HOME \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

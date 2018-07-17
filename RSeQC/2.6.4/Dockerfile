#### BASE IMAGE #####
FROM python:2.7

##### METADATA #####
LABEL base.image="python:2.7"
LABEL version="1"
LABEL software="2.6.4"
LABEL software.version="2.6.4"
LABEL software.description="An RNA-seq Quality Control Package"
LABEL software.website="http://rseqc.sourceforge.net/"
LABEL software.documentation="http://rseqc.sourceforge.net/"
LABEL software.license="http://rseqc.sourceforge.net/#license"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.6.4

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y gcc zlib1g-dev liblzo2-dev \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
  && pip install RSeQC==2.6.4

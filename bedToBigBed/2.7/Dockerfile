##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="1"
LABEL software="bedToBigBed"
LABEL software.version="2.7"
LABEL software.description="Convert bed file to bigBed. (BigBed version: 4)"
LABEL software.website="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed"
LABEL software.documentation="http://genome.ucsc.edu/goldenPath/help/bigBed.html"
LABEL software.license="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/"
LABEL software.tags="Genomics"
LABEL maintainer="christina.herrmann@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.7

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget libkrb5-3 \
  && wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed \
  && chmod +x bedToBigBed \
  && mv bedToBigBed /usr/bin \
# Clean up
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \

##### WORKING DIRECTORY #####
WORKDIR /data/

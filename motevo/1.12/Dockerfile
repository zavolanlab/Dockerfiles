##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="motevo"
LABEL software.version="1.12"
LABEL software.description="Inference of regulatory motifs from multiple alignments."
LABEL software.website="http://nimwegenlab.org/?page_id=133"
LABEL software.documentation="http://nimwegenlab.org/?page_id=133"
LABEL software.license=""
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.12
ENV DEBIAN_FRONTEND noninteractive

COPY motevo_ver${SOFTWARE_VERSION}.tar.gz $HOME

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes make g++ libgslcblas0 libglib2.0-dev libgtk2.0 \
  && tar xzvf motevo_ver${SOFTWARE_VERSION}.tar.gz \
  && cd motevo_ver${SOFTWARE_VERSION}/source \
  && make \
  && cd ../ \
  && cp bin/* /usr/local/bin/ \
  && cd ../ \
  && rm motevo_ver${SOFTWARE_VERSION}.tar.gz \
  && apt-get remove --purge --yes make g++  \
  && apt-get autoremove --purge --yes

##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="1"
LABEL software="segemehl"
LABEL software.version="0.2.0"
LABEL software.description="short read mapping with gaps"
LABEL software.website="http://www.bioinf.uni-leipzig.de/Software/segemehl/"
LABEL software.documentation="http://www.bioinf.uni-leipzig.de/Software/segemehl/"
LABEL software.license="link to license(file) of original software"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 0.2.0

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget make gcc build-essential zlib1g-dev libncurses5-dev samtools msmtp zlib1g-dev make libncurses5-dev libxml2-dev \
  && wget http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_0_2_0.tar.gz \
  && tar xzvf segemehl_0_2_0.tar.gz \
  && cd segemehl_0_2_0/segemehl/ \
  && make \
  && cd ../.. \
  && cp segemehl_0_2_0/segemehl/segemehl.x /usr/local/bin \
  && cp segemehl_0_2_0/segemehl/lack.x /usr/local/bin \
  && rm -rf segemehl_0_2_0* segemehl_0_2_0* \
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="2"
LABEL software="GMAP"
LABEL software.version="2018-05-11"
LABEL software.description="GMAP: A Genomic Mapping and Alignment Program for mRNA and EST Sequences, and GSNAP: Genomic Short-read Nucleotide Alignment Program"
LABEL software.website="http://research-pub.gene.com/gmap/"
LABEL software.documentation="http://research-pub.gene.com/gmap/"
LABEL software.license="link to license(file) of original software"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2018-05-11

##### INSTALL #####

RUN apt-get update -y \
  && apt-get install --y make g++ zlib1g-dev libbz2-dev wget perl \
  && cd /tmp \
  && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-05-11.tar.gz \
  && tar xzf gmap-gsnap-*.tar.gz && \
  && cd /tmp/gmap-* && \
  && ./configure --prefix=/usr/local --with-simd-level=sse42 \
  && make \
  && make check \
  && make install \
  && rm -rf /tmp/gmap-* \
  && apt-get remove --purge --yes make \
  && apt-get remove --purge --yes g++ \
  && apt-get remove --purge --yes wget \
  && apt-get autoremove --purge --yes

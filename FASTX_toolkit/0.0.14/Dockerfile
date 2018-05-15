##### BASE IMAGE #####
FROM ubuntu:16.04

##### METADATA #####
LABEL base.image="ubuntu:16.04"
LABEL version="1"
LABEL software="FASTX_toolkit"
LABEL software.version="0.0.14"
LABEL software.description="The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. "
LABEL software.website="http://hannonlab.cshl.edu/fastx_toolkit/"
LABEL software.documentation="http://hannonlab.cshl.edu/fastx_toolkit/commandline.html"
LABEL software.license="http://hannonlab.cshl.edu/fastx_toolkit/license.html"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="christina.herrmann@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV libgtextutils_version 0.7
ENV fastx_version 0.0.14

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install --yes gcc g++ pkg-config wget \
  && wget https://github.com/agordon/libgtextutils/releases/download/${libgtextutils_version}/libgtextutils-${libgtextutils_version}.tar.gz \
  && tar -xzf libgtextutils-${libgtextutils_version}.tar.gz \
  && cd libgtextutils-${libgtextutils_version} \
  && ./configure \
  && make \
  && make install \
  && cd .. \
  && wget https://github.com/agordon/fastx_toolkit/releases/download/${fastx_version}/fastx_toolkit-${fastx_version}.tar.bz2 \
  && tar -xjf fastx_toolkit-${fastx_version}.tar.bz2 \
  && cd fastx_toolkit-${fastx_version} \
  && ./configure \
  && make \
  && make install \
  && cd .. \
# Remove files no longer needed
  && rm -r libgtextutils-${libgtextutils_version}.tar.gz  libgtextutils-${libgtextutils_version} fastx_toolkit-${fastx_version}.tar.bz2 fastx_toolkit-${fastx_version}\
# Clean up
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \


##### WORKING DIRECTORY #####
WORKDIR /data/

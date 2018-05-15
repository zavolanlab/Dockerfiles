##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="1"
LABEL software="STAR"
LABEL software.version="2.6.0a"
LABEL software.description="Spliced Transcripts Alignment to a Reference"
LABEL software.website="https://github.com/alexdobin/STAR"
LABEL software.documentation="https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
LABEL software.license="https://github.com/alexdobin/STAR/blob/master/LICENSE"
LABEL software.tags="Transcriptomics"
LABEL maintainer="christina.herrmann@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.6.0a

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget python build-essential libz-dev\
  && wget https://github.com/alexdobin/STAR/archive/${SOFTWARE_VERSION}.tar.gz \
  && tar -xzf ${SOFTWARE_VERSION}.tar.gz \
  && cd STAR-${SOFTWARE_VERSION}/source \
  && make STAR \
  && cp STAR ../../usr/local/bin \
  && cd ../.. \
  && rm -r STAR-${SOFTWARE_VERSION} ${SOFTWARE_VERSION}.tar.gz \
# Clean up
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \


##### WORKING DIRECTORY #####
WORKDIR /data/

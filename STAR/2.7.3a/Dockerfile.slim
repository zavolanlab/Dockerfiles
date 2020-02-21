##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL software="STAR"
LABEL software.version="2.7.3a"
LABEL software.description="Spliced Transcripts Alignment to a Reference"
LABEL software.website="https://github.com/alexdobin/STAR"
LABEL software.documentation="https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
LABEL software.license="https://github.com/alexdobin/STAR/blob/master/LICENSE"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"
LABEL base.image="ubuntu:18.04"

##### VARIABLES #####
ENV SOFTWARE_VERSION 2.7.3a
ENV PACKAGES wget ca-certificates python build-essential libz-dev

##### INSTALLATION #####
RUN apt-get update
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget https://github.com/alexdobin/STAR/archive/${SOFTWARE_VERSION}.tar.gz
RUN tar -xzvf ${SOFTWARE_VERSION}.tar.gz
RUN cd STAR-${SOFTWARE_VERSION}/source && make STAR STARlong

#### CLEAN IMAGE ####
FROM ubuntu:18.04
ENV SOFTWARE_VERSION 2.7.3a
ENV PACKAGES libgomp1
COPY --from=build /STAR-${SOFTWARE_VERSION}/source/STAR* /usr/bin/
RUN apt-get update \
  && apt-get install -y --no-install-recommends ${PACKAGES} \
  && rm -rf /var/lib/apt/lists/*


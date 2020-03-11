##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL software="gffcompare"
LABEL software.version="0.11.5"
LABEL software.description="classify, merge, tracking and annotation of GFF files by comparing to a reference annotation GFF"
LABEL software.website="https://github.com/gpertea/gffcompare"
LABEL software.documentation="https://github.com/gpertea/gffcompare"
LABEL software.license="https://github.com/gpertea/gffcompare/blob/master/LICENSE"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"
LABEL base.image="ubuntu:18.04"

##### VARIABLES #####
ENV SOFTWARE_VERSION 0.11.5
ENV PACKAGES wget ca-certificates

##### INSTALLATION #####
RUN apt-get update
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-${SOFTWARE_VERSION}.Linux_x86_64.tar.gz
RUN tar -xzvf gffcompare-${SOFTWARE_VERSION}.Linux_x86_64.tar.gz
RUN ls /

#### CLEAN IMAGE ####
FROM ubuntu:18.04
ENV SOFTWARE_VERSION 0.11.5
COPY --from=build /gffcompare-${SOFTWARE_VERSION}.Linux_x86_64/gffcompare /gffcompare-${SOFTWARE_VERSION}.Linux_x86_64/trmap /usr/bin/

##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL software="gffread"
LABEL software.version="0.11.7"
LABEL software.description="GFF/GTF utility providing format conversions, region filtering, FASTA sequence extraction and more"
LABEL software.website="https://github.com/gpertea/gffread"
LABEL software.documentation="https://github.com/gpertea/gffread"
LABEL software.license="https://github.com/gpertea/gffread/blob/master/LICENSE"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"
LABEL base.image="ubuntu:18.04"

##### VARIABLES #####
ENV SOFTWARE_VERSION 0.11.7
ENV PACKAGES wget ca-certificates

##### INSTALLATION #####
RUN apt-get update
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget https://github.com/gpertea/gffread/releases/download/v0.11.7/gffread-0.11.7.Linux_x86_64.tar.gz
RUN tar -xzvf gffread-${SOFTWARE_VERSION}.Linux_x86_64.tar.gz

#### CLEAN IMAGE ####
FROM ubuntu:18.04
ENV SOFTWARE_VERSION 0.11.7
COPY --from=build /gffread-${SOFTWARE_VERSION}.Linux_x86_64/gffread /usr/bin/
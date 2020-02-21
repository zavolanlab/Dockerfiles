##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="0.46.1"
LABEL software="kallisto"
LABEL software.version="0.46.1"
LABEL software.description="Quantification of transcripts expression"
LABEL software.website="https://pachterlab.github.io/kallisto/about"
LABEL software.documentation=""
LABEL software.license=""
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

#### VARIABLES #####
ENV SOFTWARE_VERSION 0.46.1
ENV PACKAGES wget ca-certificates

#### INSTALL ####
RUN apt-get update
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget https://github.com/pachterlab/kallisto/releases/download/v${SOFTWARE_VERSION}/kallisto_linux-v${SOFTWARE_VERSION}.tar.gz
RUN tar -xzvf kallisto_linux-v${SOFTWARE_VERSION}.tar.gz

#### CREATE CLEAN IMAGE ####
FROM ubuntu:18.04
COPY --from=build /kallisto/kallisto /usr/bin/


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

RUN apt-get update \
  && apt-get install -y wget \
  && wget https://github.com/pachterlab/kallisto/releases/download/v${SOFTWARE_VERSION}/kallisto_linux-v${SOFTWARE_VERSION}.tar.gz \
  && tar -xzvf kallisto_linux-v${SOFTWARE_VERSION}.tar.gz \
  && mv /kallisto/kallisto /usr/bin/

FROM ubuntu:18.04
COPY --from=build /usr/bin/kallisto /usr/bin/kallisto


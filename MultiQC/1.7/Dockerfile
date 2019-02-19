##### BASE IMAGE #####
FROM python:3.6.5

##### METADATA #####
LABEL base.image="python:3.6.5"
LABEL version="1"
LABEL software="multiqc"
LABEL software.version="1.7"
LABEL software.description="MultiQC is a tool to aggregate bioinformatics results across many samples into a single report."
LABEL software.website="https://multiqc.info/"
LABEL software.documentation="https://multiqc.info/docs/"
LABEL software.license="GNU GPLv3"
LABEL software.tags="General,Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.7

##### INSTALL #####
RUN apt-get update -y \
  && pip install --upgrade pip \
  && pip install multiqc==${SOFTWARE_VERSION} \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \


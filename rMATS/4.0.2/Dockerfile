##### BASE IMAGE #####

FROM continuumio/miniconda:4.5.12

##### METADATA #####
LABEL base.image="continuumio/miniconda:4.5.12"
LABEL version="1"
LABEL software="rmats"
LABEL software.version="4.0.2"
LABEL software.description="MATS is a computational tool to detect differential alternative splicing events from RNA-Seq data."
LABEL software.website="http://rnaseq-mats.sourceforge.net/"
LABEL software.documentation="http://rnaseq-mats.sourceforge.net/user_guide.htm"
LABEL software.license=""
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 4.0.2

##### INSTALL #####
RUN apt-get update -y \
  && conda install -c bioconda rmats=4.0.2

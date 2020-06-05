##### BASE IMAGE #####
FROM python:3.6.6-slim

##### METADATA #####
LABEL base.image="python:3.6.6-slim"
LABEL version="1"
LABEL software="zgtf"
LABEL software.version="0.1"
LABEL software.description="gtf conversion utlity"
LABEL software.website="https://github.com/zavolanlab/zgtf"
LABEL software.documentation="https://github.com/zavolanlab/zgtf"
LABEL software.license="https://github.com/zavolanlab/zgtf/blob/master/LICENSE"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="Apache License 2.0"

##### INSTALL #####
RUN apt-get update \
  && git clone -branch 0.1 https://github.com/zavolanlab/zgtf.git \
  && pip install .

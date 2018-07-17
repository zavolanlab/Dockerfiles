##### BASE IMAGE #####
FROM python:2.7.15

##### METADATA #####
LABEL base.image="python:2.7.15"
LABEL version="1"
LABEL software="2.7.15"
LABEL software.version="2.7.15"
LABEL software.description="Python programming language"
LABEL software.website="https://www.python.org/"
LABEL software.documentation="https://docs.python.org/3/"
LABEL software.license="https://docs.python.org/3/license.html"
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.7.15

##### INSTALL #####
RUN apt-get update -y \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

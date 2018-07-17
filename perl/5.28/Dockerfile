##### BASE IMAGE #####
FROM perl:5.28

##### METADATA #####
LABEL base.image="perl.5.28"
LABEL version="1"
LABEL software="5.28"
LABEL software.version="5.28"
LABEL software.description="Perl is a high-level, general-purpose, interpreted, dynamic programming language."
LABEL software.website="https://www.perl.org/"
LABEL software.documentation="https://www.perl.org/docs.html"
LABEL software.license=""
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 5.28

##### INSTALL #####
RUN apt-get update -y \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

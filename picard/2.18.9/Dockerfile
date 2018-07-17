##### BASE IMAGE #####
FROM java:8

##### METADATA #####
LABEL base.image="java:8"
LABEL version="1"
LABEL software="2.18.9"
LABEL software.version="2.18.9"
LABEL software.description="Command line tools for manipulating high-throughput sequencing (HTS) data and formats"
LABEL software.website="http://broadinstitute.github.io/picard/"
LABEL software.documentation="http://broadinstitute.github.io/picard/command-line-overview.html"
LABEL software.license="https://github.com/broadinstitute/picard/blob/master/LICENSE.txt"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.18.9

##### INSTALL #####
RUN apt-get update \
    && apt-get --no-install-recommends install -y --force-yes git ant wget \
    && apt-get clean autoclean \
    && apt-get autoremove -y \
    && wget https://github.com/broadinstitute/picard/releases/download/${SOFTWARE_VERSION}/picard.jar \
    && cp picard.jar /usr/local/bin/

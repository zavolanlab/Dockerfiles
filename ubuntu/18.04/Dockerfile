##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="ubuntu"
LABEL software.version="18.04"
LABEL software.description="Ubuntu 18.04 with some basic packages (e.g. wget)"
LABEL software.website="https://www.ubuntu.com/"
LABEL software.documentation="https://www.ubuntu.com/"
LABEL software.license="https://www.ubuntu.com/licensing"
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

RUN apt-get update \
  && apt-get install -y wget

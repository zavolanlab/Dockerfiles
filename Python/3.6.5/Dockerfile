##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="1"
LABEL software="Python3"
LABEL software.version="3.6.5"
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
ENV SOFTWARE_VERSION 3.6.5

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget build-essential libz-dev libssl-dev openssl gcc \
  && wget https://www.python.org/ftp/python/3.6.5/Python-3.6.5.tgz \
  && tar zxf Python-3.6.5.tgz \
  && cd Python-3.6.5 \
  && ./configure --prefix=/usr/local \
  && make \
  && make altinstall \
  && ln -s /usr/local/bin/python3.6 /usr/local/bin/python \
  && ln -s /usr/local/bin/pip3.6 /usr/local/bin/pip \
  && cd .. \
  && rm -r Python-3.6.5.tgz \
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

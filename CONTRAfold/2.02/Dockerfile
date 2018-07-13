FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="CONTRAfold"
LABEL software.version="2.02"
LABEL software.description="CONTRAfold is a novel secondary structure prediction method based on conditional log-linear models (CLLMs)."
LABEL software.website="http://contra.stanford.edu/contrafold"
LABEL software.documentation="http://contra.stanford.edu/contrafold/manual_v2_02.pdf"
LABEL software.license="SPDX:BSD-2-Clause"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="jaroslaw.surkont@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

RUN apt-get update && \
  apt-get install --yes \
    make \
    g++ \
    wget && \
  cd /tmp && \
  wget http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz && \
  tar xzf contrafold_*.tar.gz && \
  cd /tmp/contrafold/src && \
  sed -i '/^CXXFLAGS = / s/$/ -fpermissive/' Makefile && \
  sed -i '/#include <vector>/a #include <limits.h>' Utilities.hpp && \
  make && \
  mv contrafold /usr/bin/ && \
  rm -rf /tmp/contrafold* && \
  apt-get remove --purge --yes \
    make \
    g++ \
    wget && \
  apt-get autoremove --purge --yes

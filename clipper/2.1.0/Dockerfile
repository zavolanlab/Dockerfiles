##### BASE IMAGE #####
FROM python:3.8-slim

##### METADATA #####
LABEL base.image="python:3.8-slim"
LABEL version="1"
LABEL software="clipper"
LABEL software.version="2.1.0"
LABEL software.description="CLIPper - CLIP peak enrichment recognition"
LABEL software.website="https://github.com/YeoLab/clipper"
LABEL software.documentation="https://github.com/YeoLab/clipper/wiki/CLIPper-Home"
LABEL software.license="GPL2"
LABEL software.tags="Transcriptomics, Clip, RBPs"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.1.0

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y git gcc g++ zlib1g-dev bedtools \
  && git clone https://github.com/YeoLab/clipper.git \
  && cd clipper \
  && git checkout 8566530fa8b27b8840479fde6860572007c88261 \
  && pip install numpy==1.18.5 \
  && pip install cython pysam>=0.15.3 scipy>=1.5.0 matplotlib>=3.2.2 pybedtools scikit-learn>=0.23.1 HTSeq>=0.11.3 \
  && pip install . \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

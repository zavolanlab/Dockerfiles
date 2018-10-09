##### BASE IMAGE #####
FROM python:2.7.15-slim

##### METADATA #####
LABEL base.image="python:2.7.15-slim"
LABEL version="1"
LABEL software="clipper"
LABEL software.version="0.2.0"
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
ENV SOFTWARE_VERSION 0.2.0

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y git \
  && git clone https://github.com/YeoLab/clipper.git \
  && cd clipper \
  && git checkout 423b54b9995d248197a153a066c880348ebd4eef

RUN apt-get install -y gcc g++ libxml2-dev libxslt1-dev zlib1g-dev libpng-dev libfreetype6-dev \
  && pip install numpy==1.5.1 \
  && pip install cython pandas==0.22.0 matplotlib==1.1.0  pysam>=0.9.0 scipy pybedtools scikit-learn==0.17.0 HTSeq \
  && cd clipper \
  && python setup.py install \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



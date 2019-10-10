###############################################################################
#
#   Dockerfile for the container
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 06-10-2019
#   LICENSE: GPL
#
###############################################################################

# BASE IMAGE
FROM ubuntu:18.04

# PATH EXPORT
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

# METADATA
LABEL base.image="ubuntu:18.04"
LABEL version="0.0.9"
LABEL software="MAPP"
LABEL software.version="0.0.9"
LABEL software.description=\
"Dependencies for the pipeline to analyse motif activity on alternative \
splicing and polyadenylation."
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation="Swiss Institute of Bioinformatics"
LABEL maintainer.location=\
"Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"

# UPDATE THE PACKAGE LIST
RUN apt-get update

# INSTALL WGET
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# INSTALL MINICONDA WITH PYTHON 3.7
RUN wget --no-verbose \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

# UPDATE CONDA
RUN conda update --name base --channel defaults conda

# COPY THE YAML & INSTALL SOFTWARE WITH CONDA
COPY packages.yaml .
RUN conda env update --name base --file packages.yaml
RUN conda clean -all

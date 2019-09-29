##### BASE #####
FROM rocker/r-ver:3.6.0

##### METADATA #####
MAINTAINER alexander.kanitz@alumni.ethz.ch
LABEL base.image="rocker/r-ver:3.6.0"
LABEL version="2"
LABEL software="ascii_alignment_pileup.R"
LABEL software.version="1.0.1"
LABEL software.description="Generates an ASCII-style pileup of read alignments in one or more BAM files against one or more regions specified in a BED file."
LABEL software.website="https://github.com/uniqueg/scripts/tree/master/scripts/ascii_alignment_pileup"
LABEL software.documentation="https://github.com/uniqueg/scripts/tree/master/scripts/ascii_alignment_pileup/README.md"
LABEL software.license="https://spdx.org/licenses/Apache-2.0"
LABEL software.tags="BED,BAM,GFF,FASTA,alignments,pileup,ASCII,microRNA,miRNA,miR,isomiR"
LABEL maintainer="alexander.kanitz@alumni.ethz.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"

##### INSTALL #####
RUN \
  apt-get update -y \
  && apt-get install -y --no-install-recommends \
    apt-transport-https \
    curl \
    libbamtools-dev \
    libbz2-dev \
    libcurl4-gnutls-dev \
    liblzma-dev \
    libssl-dev \
    libxml2-dev \
    software-properties-common \
    zlib1g-dev \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN \
  install2.r \
    --error \
    --skipinstalled \
    -r "https://cran.rstudio.com" \
    -r "http://www.bioconductor.org/packages/release/bioc" \
    -r "https://bioconductor.org/packages/3.9/data/annotation" \
    optparse \
    rtracklayer \
  && rm -rf /tmp/downloaded_packages/

##### VARIABLES #####
ARG WORKDIR="/home/bioc/"
ARG USER="bioc"
ARG GROUP="bioc"
ENV PATH="${WORKDIR}:${PATH}"

##### COPY #####
COPY ./ascii_alignment_pileup.R ${WORKDIR}
COPY ./tests ${WORKDIR}/tests

##### USER #####
RUN groupadd -r ${GROUP} && useradd --no-log-init -r -g ${GROUP} ${USER}

##### ENVIRONMENT  #####
WORKDIR ${WORKDIR}
RUN chown -R ${USER}:${GROUP} ${WORKDIR} \
  && chmod 700 ${WORKDIR}
USER ${USER}


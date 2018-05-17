##### BASE IMAGE #####
FROM debian:9.3

##### METADATA #####
LABEL base.image="debian:9.3"
LABEL version="1"
LABEL software="cufflinks"
LABEL software.version="2.2.1"
LABEL software.description="Transcriptome assembly and differential expression analysis for RNA-Seq."
LABEL software.website="http://cole-trapnell-lab.github.io/cufflinks/"
LABEL software.documentation="http://cole-trapnell-lab.github.io/cufflinks/manual/"
LABEL software.license="OSI-approved Boost License"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 2.2.1

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget python build-essential libz-dev \
  && wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/links-2.2.1.Linux_x86_64.tar.gz \
  && tar xzvf cufflinks-2.2.1.Linux_x86_64.tar.gz \
  && cd cufflinks-2.2.1.Linux_x86_64 \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cuffcompare /usr/local/bin/cuffcompare \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cuffdiff /usr/local/bin/cuffdiff \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cufflinks /usr/local/bin/cufflinks \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cuffmerge /usr/local/bin/cuffmerge \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cuffnorm /usr/local/bin/cuffnorm \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/cuffquant /usr/local/bin/cuffquant \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/gffread /usr/local/bin/gffread \
  && ln -s /cufflinks-2.2.1.Linux_x86_64/gtf_to_sam /usr/local/bin/gtf_to_sam \
  && cd ../ \
  && rm -rf cufflinks-2.2.1.Linux_x86_64.tar.gz \
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

##### BASE IMAGE #####
FROM python:3.6.5

##### METADATA #####
LABEL base.image="python:3.6.5"
LABEL version="2"
LABEL software="3.6.5_0.10.0_1.7.1"
LABEL software.version="3.6.5_0.10.0_1.7.1"
LABEL software.description="Python programming language|Analysing high-throughput sequencing data with Python|Biopython is a set of freely available tools for biological computation written in Python by an international team of developers."
LABEL software.website="https://www.python.org/|https://htseq.readthedocs.io/en/release_0.10.0/|https://biopython.org/"
LABEL software.documentation="https://docs.python.org/3/|https://htseq.readthedocs.io/en/release_0.10.0|https://biopython.org/wiki/Documentation"
LABEL software.license="https://docs.python.org/3/license.html|https://htseq.readthedocs.io/en/release_0.10.0/overview.html#license|https://github.com/biopython/biopython/blob/master/LICENSE.rst"
LABEL software.tags="General|Genomics,Transcriptomics|Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.6.5_0.10.0_1.7.1

##### INSTALL #####
RUN apt-get update -y \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
  && pip install HTSeq==0.10.0 \
  && pip install biopython==1.71


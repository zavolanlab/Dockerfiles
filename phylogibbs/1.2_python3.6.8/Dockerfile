##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="phylogibbs"
LABEL software.version="1.2"
LABEL software.description="Discover regulatory sites in a collection of DNA sequences."
LABEL software.website="http://www.phylogibbs.unibas.ch/cgi-bin/phylogibbs.pl"
LABEL software.documentation="http://www.phylogibbs.unibas.ch/cgi-bin/phylogibbs.pl?part=documentation&page=phylogibbs"
LABEL software.license=""
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.2
ENV LIBRARY_VERSION 1.16
ENV DEBIAN_FRONTEND noninteractive

COPY phylogibbs-${SOFTWARE_VERSION}.zip $HOME
COPY gsl-${LIBRARY_VERSION}.tar.gz $HOME

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes make g++ unzip wget libgslcblas0 libglib2.0-dev libgtk2.0 curl \
  && tar xzvf gsl-${LIBRARY_VERSION}.tar.gz \
  && cd gsl-${LIBRARY_VERSION} \
  && ./configure \
  && make \
  && make install \
  && cd ../ \
  && unzip phylogibbs-${SOFTWARE_VERSION}.zip \
  && cd phylogibbs-${SOFTWARE_VERSION} \
  && ./configure \
  && make \
  && make install \
  && cd ../ \
  && rm phylogibbs-${SOFTWARE_VERSION}.zip \
  && rm gsl-${LIBRARY_VERSION}.tar.gz \
  && unlink /usr/bin/python \
  && ln -s /usr/bin/python3 /usr/bin/python \
  && curl https://bootstrap.pypa.io/ez_setup.py -o - | python \
  && rm setuptools-* \
  && pip install pandas==0.23.3 seaborn==0.9.0 pysam==0.14.1 pyfasta==0.5.2 scikit-learn==0.19.2 \
  && apt-get remove --purge --yes make g++ unzip \
  && apt-get autoremove --purge --yes

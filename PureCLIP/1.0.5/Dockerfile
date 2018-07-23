##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="PureClip"
LABEL software.version="1.0.5"
LABEL software.description="PureCLIP: capturing target-specific protein–RNA interaction footprints from single-nucleotide CLIP-seq data"
LABEL software.website="https://github.com/skrakau/PureCLIP"
LABEL software.documentation="https://github.com/skrakau/PureCLIP"
LABEL software.license="https://github.com/skrakau/PureCLIP/blob/master/LICENSE.md"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"


##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.0.5
ENV LIBRARY_VERSION 1.16
ENV DEBIAN_FRONTEND noninteractive

COPY v${SOFTWARE_VERSION}.tar.gz $home
COPY gsl-${LIBRARY_VERSION}.tar.gz $HOME

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes make g++ cmake libgslcblas0 libglib2.0-dev libgtk2.0 \
  && tar xzvf gsl-${LIBRARY_VERSION}.tar.gz \
  && cd gsl-${LIBRARY_VERSION} \
  && ./configure \
  && make \
  && make install \
  && cd ../ \
  && tar xzvf v${SOFTWARE_VERSION}.tar.gz \
  && cd PureCLIP-${SOFTWARE_VERSION} \
  && mkdir build \
  && cd build \
  && cmake ../src \
  && make \
  && cp pureclip /usr/local/bin/ \
  && cp winextract /usr/local/bin/ \ 
  && cd ../../ \
  && rm v${SOFTWARE_VERSION}.tar.gz \
  && rm gsl-${LIBRARY_VERSION}.tar.gz \
  && apt-get remove --purge --yes make g++  \
  && apt-get autoremove --purge --yes

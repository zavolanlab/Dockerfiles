##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="omniCLIP"
LABEL software.version="master_77159cd1a0119cbf3b067dbfed0cba73b1091e91"
LABEL software.description="omniCLIP is a CLIP-Seq peak caller"
LABEL software.website="https://github.com/philippdre/omniCLIP"
LABEL software.documentation="https://github.com/philippdre/omniCLIP"
LABEL software.license="https://github.com/philippdre/omniCLIP/blob/master/LICENSE"
LABEL software.tags="Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

WORKDIR /opt/

##### INSTALL #####
RUN apt-get update && \
  apt-get install -y tzdata && \
  ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime && \
  dpkg-reconfigure --frontend noninteractive tzdata && \
  apt-get install -y git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev ca-certificates zlib1g-dev curl unzip autoconf python-pip python-numpy python-tk && \
  pip install \
    biopython==1.68 \
    brewer2mpl==1.4 \
    cython==0.24.1 \
    gffutils==0.8.7.1 \
    h5py==2.6.0 \
    intervaltree==2.1.0 \
    matplotlib==1.5.3 \
    pandas==0.19.2 \
    prettyplotlib==0.1.7 \
    pysam==0.9.1.4 \
    scikit-learn==0.18.1 \
    scipy==0.19.0 \
    statsmodels==0.6.1

RUN cd /opt/ && \
  git clone https://github.com/philippdre/omniCLIP.git && \
  cd omniCLIP && \
  git checkout 77159cd1a0119cbf3b067dbfed0cba73b1091e91 && \
  chmod +x /opt/omniCLIP/omniCLIP.py && \
  chmod +x /opt/omniCLIP/data_parsing/CreateGeneAnnotDB.py && \
  rm -rf .git && \
  cd /opt/omniCLIP/stat/ && \
  ./CompileCython.sh

#  && apt-get autoremove -y \
#  && apt-get clean \
#  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /opt/omniCLIP/omniCLIP.py /usr/local/bin/omniCLIP.py && \
    ln -s /opt/omniCLIP/data_parsing/CreateGeneAnnotDB.py /usr/local/bin/CreateGeneAnnotDB.py

ENV PYTHONPATH=${PYTHONPATH}:/opt/omniCLIP:/opt/omniCLIP/stat:/opt/omniCLIP/data_parsing

WORKDIR /opt/omniCLIP/

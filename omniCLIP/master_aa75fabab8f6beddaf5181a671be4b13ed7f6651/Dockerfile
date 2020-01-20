##### BASE IMAGE #####
FROM python:2.7.17

##### METADATA #####
LABEL base.image="python:2.7.17"
LABEL version="1"
LABEL software="omniCLIP"
LABEL software.version="master_aa75fabab8f6beddaf5181a671be4b13ed7f6651"
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
  apt-get install -y python-numpy && \  
  pip install \
    biopython==1.70 \
    brewer2mpl==1.4.1 \
    Cython==0.28.1 \
    gffutils==0.9 \
    h5py==2.7.1 \
    intervaltree==2.1.0 \
    matplotlib==2.2.2 \
    multiprocessing==2.6.2.1 \
    pandas==0.22.0 \
    prettyplotlib==0.1.7 \
    pysam==0.14 \
    scikit-learn==0.19.1 \
    scipy==1.0.0 \
    statsmodels==0.8.0

RUN cd /opt/ && \
  git clone https://github.com/philippdre/omniCLIP.git && \
  cd omniCLIP && \
  git checkout aa75fabab8f6beddaf5181a671be4b13ed7f6651 && \
  chmod +x /opt/omniCLIP/omniCLIP.py && \
  chmod +x /opt/omniCLIP/data_parsing/CreateGeneAnnotDB.py && \
  rm -rf .git && \
  cd /opt/omniCLIP/stat/ && \
  sed -i 's/\/usr\/include\/python2.7/\/usr\/local\/include\/python2.7/' CompileCython.sh && \
  ./CompileCython.sh

ENV PYTHONPATH=${PYTHONPATH}:/opt/omniCLIP:/opt/omniCLIP/stat:/opt/omniCLIP/data_parsing

WORKDIR /opt/omniCLIP/

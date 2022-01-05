##### BASE IMAGE #####
FROM python:3.8

##### METADATA #####
LABEL base.image="python:3.8"
LABEL version="1"
LABEL software="omniCLIP"
LABEL software.version="master_83e44f8a044227afb3d457ee33a8da7e07964d91_5"
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

##### INSTALL #####
RUN apt-get update && \
  apt-get install git && \
  pip install Cython numpy

RUN git clone -b omniclip_fixes https://github.com/fgypas/omniCLIP.git omniCLIP_clone && \
  cd omniCLIP_clone && \
  git checkout 83e44f8a044227afb3d457ee33a8da7e07964d91 && \
  rm -rf .git && \
  cd ../ && \
  mv -v /omniCLIP_clone/* . && \
  rm -rf omniCLIP_clone/ && \
  chmod +x /omniCLIP/omniCLIP.py && \
  python3 setup.py install
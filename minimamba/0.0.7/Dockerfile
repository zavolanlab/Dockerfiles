##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="minimamba"
LABEL software.version="0.0.7"
LABEL software.description="Mamba is a fast cross platform package manager, based on and compatible with conda."
LABEL software.website="https://quantstack.net/mamba.html"
LABEL software.documentation="https://github.com/QuantStack/mamba"
LABEL software.license="https://github.com/QuantStack/mamba/blob/master/LICENSE"
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes wget \
  && wget https://github.com/QuantStack/mamba/releases/download/0.0.7/minimamba-0.0.7-Linux-x86_64-py37.sh -O ~/minimamba.sh \
  && bash ~/minimamba.sh -b -p /opt/conda \
  && rm ~/minimamba.sh \
  && apt-get remove --purge --yes wget \
  && apt-get autoremove --purge --yes

ENV PATH /opt/conda/bin:${PATH}

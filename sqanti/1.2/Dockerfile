##### BASE IMAGE #####
FROM ubuntu:16.04

##### METADATA #####
LABEL base.image="ubuntu:16.04"
LABEL version=""
LABEL software=""
LABEL software.version=""
LABEL software.description=""
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license=""
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
#ENV SOFTWARE_VERSION 2018-05-11

##### INSTALL #####

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install -y vim wget curl zlib1g-dev git software-properties-common libcurl4-gnutls-dev libxml2-dev libssl-dev apt-transport-https libmariadb-client-lgpl-dev g++ libbz2-dev perl libgit2-dev \
  && apt-get install -y python-pip python-dev build-essential \
  && pip install Pysam==0.15 \
  && pip install psutil==5.4.7

RUN cd /tmp \
  && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-05-11.tar.gz \
  && tar xzf gmap-gsnap-*.tar.gz \
  && cd /tmp/gmap-* \
  && ./configure --prefix=/usr/local --with-simd-level=sse42 \
  && make \
  && make check \
  && make install \
  && rm -rf /tmp/gmap-* \
  && cd $HOME

#RUN apt-get install -y libgit2-dev

RUN git clone https://fgypas@bitbucket.org/ConesaLab/sqanti.git \
  && cd sqanti \
  && cp sqanti_qc.py /usr/local/bin/ \
  && cp qanti_filter.py /usr/local/bin/ \
  && cp -r utilities /usr/local/bin/

FROM r-base:3.5.1

RUN Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("ggplot2", "scales", "reshape", "gridExtra", "grid", "rpart", "ROCR", "caret", "lattice", "foreach", "e1071", "randomForest", "partykit", "ipred", "rpart.plot", "doMC", "nnet", "ROSE", "pROC", "MLmetrics"), ask=FALSE)'



#  && apt-get remove --purge --yes make g++ wget \
#  && apt-get autoremove --purge --yes

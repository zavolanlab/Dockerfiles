FROM rocker/tidyverse:3.5.3
MAINTAINER Alex Kanitz "alexander.kanitz@alumni.ethz.ch"

## Install system dependencies
RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  curl \
  git \
  libarmadillo7 \
  libarmadillo-dev \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libhdf5-cpp-100 \
  libpng-dev \
  libssl-dev \
  libxml2-dev \
  libz-dev

## Install Python
RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tar.xz \
  && tar xJf Python-3.6.0.tar.xz \
  && cd Python-3.6.0 \
  && ./configure \
  && make altinstall \
  && ln -s /Python-3.6.0/python /usr/local/bin \
  && cd .. \
  && python -m pip install --upgrade pip setuptools wheel virtualenv

## Install Python dependencies
RUN pip3 install umap-learn

## Install R packages
RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("Rcpp","RcppArmadillo", "Matrix", "mgcv", "abind","igraph","h5","Rtsne", "cluster","data.table"))'

## Install Bioconductor packages
RUN \
  R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("pcaMethods","edgeR","Rsamtools","GenomicAlignments","GenomeInfoDb","Biostrings"),suppressAutoUpdate=TRUE,ask=FALSE,suppressUpdates=TRUE)'

## Install R packages from GitHub
RUN \
  R -e 'library("devtools"); install_github(c("satijalab/seurat@develop", "satijalab/seurat-wrappers"))'

## Set user
RUN useradd -m user
USER user
ENTRYPOINT ["/bin/bash"]
WORKDIR "/home/user"

## Install Velocyto.R
RUN \
  git clone https://github.com/velocyto-team/velocyto.R && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4

RUN \
  echo '.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths()))' > .Rprofile && \
  R -e 'devtools::install_local("~/velocyto.R/",dep=T,upgrade_dependencies=F)'

ENV  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH \
  R_PROFILE=~/.Rprofile

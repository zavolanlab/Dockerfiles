##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="4"
LABEL software="sqanti"
LABEL software.version="1.2"
LABEL software.description="SQANTI: Structural and Quality Annotation of Novel Transcript Isoforms"
LABEL software.website="https://bitbucket.org/ConesaLab/sqanti"
LABEL software.documentation="https://bitbucket.org/ConesaLab/sqanti"
LABEL software.license="https://bitbucket.org/ConesaLab/sqanti/src/6927e53e56d238cb94157ec06788d1e4df198fe5/LICENSE"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability

##### INSTALL #####

# install dependencies and python
RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install -y vim wget curl zlib1g-dev git software-properties-common libcurl4-gnutls-dev libxml2-dev apt-transport-https libmariadb-client-lgpl-dev g++ libbz2-dev perl libgit2-dev libssl1.1 libssl-dev libgfortran3 \
  && apt-get install -y python-pip python-dev build-essential \
  && pip install Pysam==0.15 \
  && pip install psutil==5.4.7

# isntall gmap/gsnap
RUN cd /tmp \
  && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz \
  && tar xzf gmap-gsnap-*.tar.gz \
  && cd /tmp/gmap-* \
  && ./configure --prefix=/usr/local --with-simd-level=sse42 \
  && make \
  && make check \
  && make install \
  && rm -rf /tmp/gmap-* \
  && cd $HOME

# install R dependencies
RUN Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("ggplot2", "scales", "reshape", "gridExtra", "grid", "rpart", "ROCR", "caret", "lattice", "foreach", "e1071", "randomForest", "partykit", "ipred", "rpart.plot", "doMC", "nnet", "ROSE", "pROC", "MLmetrics"), ask=FALSE)'

# install sqanti
RUN git clone https://fgypas@bitbucket.org/ConesaLab/sqanti.git \
  && cd sqanti \
  && git checkout 6927e53e56d238cb94157ec06788d1e4df198fe5 \
  && chmod +x sqanti_qc.py sqanti_filter.py utilities/* \
  && cp -R ./ /usr/bin/

# quick fix for missing library
RUN cd /usr/lib/x86_64-linux-gnu/ \
  && ln -s libpng16.so libpng12.so.0 \
  && cd $HOME

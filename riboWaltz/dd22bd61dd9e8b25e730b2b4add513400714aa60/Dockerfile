##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="1"
LABEL software="riboWaltz"
LABEL software.version="dd22bd61dd9e8b25e730b2b4add513400714aa60"
LABEL software.description="riboWaltz: Optimization of ribosome P-site positioning in ribosome profiling data"
LABEL software.website="https://github.com/LabTranslationalArchitectomics/riboWaltz"
LABEL software.documentation="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/ReferenceManual.pdf"
LABEL software.license="https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/master/LICENSE"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### INSTALL #####

# install dependencies
RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install -y vim wget curl zlib1g-dev git software-properties-common libcurl4-gnutls-dev libxml2-dev apt-transport-https libmariadb-client-lgpl-dev g++ libbz2-dev perl libgit2-dev libssl1.1 libssl-dev libgfortran3 libopenblas-base

# install R dependencies
RUN Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'library(devtools); install_github("LabTranslationalArchitectomics/riboWaltz", dependencies=TRUE, ref="dd22bd61dd9e8b25e730b2b4add513400714aa60")'

# quick fix for missing library
RUN cd /usr/lib/x86_64-linux-gnu/ \
  && ln -s libpng16.so libpng12.so.0 \
  && cd $HOME

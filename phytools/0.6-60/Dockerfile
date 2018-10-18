##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="1"
LABEL software="phytools"
LABEL software.version="0.6-60"
LABEL software.description="https://cran.r-project.org/web/packages/phytools/"
LABEL software.website="https://cran.r-project.org/web/packages/phytools/"
LABEL software.documentation="https://cran.r-project.org/web/packages/phytools/"
LABEL software.license="https://cran.r-project.org/web/licenses/GPL-3"
LABEL software.tags="Genomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y curl zlib1g-dev software-properties-common libcurl4-gnutls-dev libxml2-dev apt-transport-https wget libssl1.1=1.1.1-1 libssl-dev=1.1.1-1 \
  && Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' \
  && Rscript -e 'install.packages(c("ape", "maps", "animation", "clusterGeneration", "coda", "combinat", "expm", "mnormt", "numDeriv", "phangorn", "plotrix", "scatterplot3d"), repos = "http://cran.us.r-project.org")' \
  && wget https://cran.r-project.org/src/contrib/phytools_0.6-60.tar.gz \
  && R CMD INSTALL phytools_0.6-60.tar.gz \
  && rm phytools_0.6-60.tar.gz

##### BASE IMAGE #####
FROM r-base:3.5.1

##### METADATA #####
LABEL base.image="r-base:3.5.1"
LABEL version="1"
LABEL software="seqlogo"
LABEL software.version="1.46.0"
LABEL software.description="seqLogo takes the position weight matrix of a DNA sequence motif and plots the corresponding sequence logo as introduced by Schneider and Stephens (1990)."
LABEL software.website="https://bioconductor.org/packages/release/bioc/html/seqLogo.html"
LABEL software.documentation="https://bioconductor.org/packages/release/bioc/vignettes/seqLogo/inst/doc/seqLogo.pdf"
LABEL software.license="LGPL (>= 2)"
LABEL software.tags="Genomics"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y curl zlib1g-dev software-properties-common libcurl4-gnutls-dev libxml2-dev apt-transport-https wget libssl1.1=1.1.1-1 libssl-dev=1.1.1-1 \
  && Rscript -e 'install.packages(c("devtools", "optparse"), repos = "http://cran.us.r-project.org")' \
  && wget https://bioconductor.org/packages/release/bioc/src/contrib/seqLogo_1.46.0.tar.gz \
  && R CMD INSTALL seqLogo_1.46.0.tar.gz \
  && rm seqLogo_1.46.0.tar.gz

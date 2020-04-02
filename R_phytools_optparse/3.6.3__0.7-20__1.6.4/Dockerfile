##### BASE IMAGE #####
FROM r-base:3.6.3

##### METADATA #####
LABEL base.image="r-base:3.6.3"
LABEL version="1"
LABEL software="R with phytools and optparse"
LABEL software.version="3.6.3|0.7-20|1.6.4"
LABEL software.description="R with phytools and optparse"
LABEL software.website="https://cran.r-project.org/web/packages/phytools/index.html|https://cran.r-project.org/web/packages/optparse/index.html"
LABEL software.documentation="https://cran.r-project.org/web/packages/phytools/index.html|https://cran.r-project.org/web/packages/optparse/index.html"
LABEL software.license="GPL-2|GPL-3"
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y build-essential gcc curl libxml2-dev libssl-dev r-cran-devtools libmagick++-dev

RUN Rscript -e 'require("devtools"); install_version("phytools", version ="0.7-20", repos = "http://cran.us.r-project.org")'

RUN Rscript -e 'require("devtools"); install_version("optparse", version ="1.6.4", repos = "http://cran.us.r-project.org")'

##### BASE IMAGE #####
FROM rocker/tidyverse:3.5.3

##### METADATA #####
LABEL base.image="rocker/tidyverse:3.5.3"
LABEL version="3"
LABEL software="3.5.3"
LABEL software.version="3.5.3"
LABEL software.description="R with tidyverse"
LABEL software.website="https://www.r-project.org/"
LABEL software.documentation="https://www.r-project.org/"
LABEL software.license=""
LABEL software.tags="General"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.5.3

##### INSTALL #####
RUN Rscript -e 'install.packages("optparse", repos = "http://cran.us.r-project.org")'

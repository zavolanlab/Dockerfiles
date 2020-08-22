##### BASE IMAGE #####
FROM python:3.7.4

##### METADATA #####
LABEL base.image="python:3.7.4"
LABEL version="2"
LABEL software="logomaker"
LABEL software.version="3.7.4"
LABEL software.description="Python libraries to plot sequence logos"
LABEL software.website="https://github.com/jbkinney/logomaker"
LABEL software.documentation="https://logomaker.readthedocs.io/en/latest/"
LABEL software.license="MIT"
LABEL software.tags="General"
LABEL maintainer="akrish136@gmail.com"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 3.7.4

##### INSTALL #####
RUN apt-get update -y \
  && pip install --upgrade pip \
  && pip install logomaker==0.8 matplotlib==3.3.0 \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
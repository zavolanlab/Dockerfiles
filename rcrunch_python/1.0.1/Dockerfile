##### BASE IMAGE #####
FROM python:3.6.6

##### METADATA #####
LABEL base.image="python:3.6.6"
LABEL version="1"
LABEL software="RCRUNCH"
LABEL software.version="1.0.1"
LABEL software.description="Python packages used by the RCRUNCH pipeline"
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license=""
LABEL software.tags="Transcriptomics"
LABEL maintainer="maria.katsantoni@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license=""

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 1.0.1

##### INSTALL #####
RUN apt-get update \
  && apt-get install -y gcc \
  && pip install pandas==0.23.3 seaborn==0.9.0 pysam==0.14.1 pyfasta==0.5.2 scikit-learn==0.19.2 ushuffle==1.1.0 Cython==0.27.1

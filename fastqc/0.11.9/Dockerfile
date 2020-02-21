##### BASE IMAGE #####
FROM java:8

##### METADATA #####
LABEL base.image="java:8"
LABEL version="1"
LABEL software="FastQC"
LABEL software.version="0.11.9"
LABEL software.description="A quality control application for high throughput sequence data"
LABEL software.website="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
LABEL software.documentation="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
LABEL software.license="GPL v3 or later."
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 0.11.9

##### INSTALL #####
RUN apt-get update \
    && apt-get --no-install-recommends install -y --force-yes wget \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${SOFTWARE_VERSION}.zip \
    && unzip fastqc_v${SOFTWARE_VERSION}.zip \
    && chmod +x /FastQC/fastqc \
    && ln -s /FastQC/fastqc /usr/local/bin/ \
    && rm fastqc_v${SOFTWARE_VERSION}.zip \
    && apt-get clean autoclean \
    && apt-get autoremove -y wget 

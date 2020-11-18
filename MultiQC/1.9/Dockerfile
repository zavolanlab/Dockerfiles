##### BASE IMAGE #####
FROM python:3.6.6

##### METADATA #####
LABEL base.image="python:3.6.6"
LABEL version="1"
LABEL software="MultiQC"
LABEL software.version="1.9"
LABEL software.description="Aggregate bioinformatics results in single report"
LABEL software.website="https://github.com/ewels/MultiQC"
LABEL software.documentation="https://multiqc.info/docs/"
LABEL software.license="https://github.com/ewels/MultiQC/blob/master/LICENSE"
LABEL software.tags="Bioinformatcs"
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"

##### INSTALL #####
RUN apt-get update \
  && pip install multiqc==1.9


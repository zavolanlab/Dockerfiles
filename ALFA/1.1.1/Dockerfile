##### BASE IMAGE #####
FROM python:3.6.6

##### METADATA #####
LABEL base.image="python:3.6.6"
LABEL version="1"
LABEL software="ALFA"
LABEL software.version="1.1.1"
LABEL software.description="ALFA provides a global overview of features distribution composing NGS dataset(s)"
LABEL software.website="https://github.com/biocompibens/ALFA"
LABEL software.documentation="https://github.com/biocompibens/ALFA#manual"
LABEL software.license="https://github.com/biocompibens/ALFA/blob/master/LICENSE.txt"
LABEL software.tags="bioinformatcs ngs"
LABEL maintainer="alexander.kanitz@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"

##### INSTALL #####
RUN apt-get update -y \
  && pip install alfa==1.1.1


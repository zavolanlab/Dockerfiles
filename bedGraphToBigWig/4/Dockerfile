##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL software="bedGraphToBigWig"
LABEL software.version="4"
LABEL software.description="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/"
LABEL software.website="https://github.com/ucscGenomeBrowser/kent"
LABEL software.documentation="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/"
LABEL software.license="custom"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="alexander.kanitz@alumni.ethz.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://github.com/ucscGenomeBrowser/kent/blob/master/LICENSE"

##### VARIABLES #####
ENV PACKAGES wget

##### INSTALLATION #####
RUN apt-get update
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig

#### CLEAN IMAGE ####
FROM ubuntu:18.04
ENV PACKAGES libkrb5-3
COPY --from=build bedGraphToBigWig /usr/bin/
RUN chmod +x /usr/bin/bedGraphToBigWig
RUN apt-get update \
  && apt-get install -y --no-install-recommends ${PACKAGES} \
  && rm -rf /var/lib/apt/lists/*


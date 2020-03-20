##### BASE IMAGE #####
FROM ubuntu:18.04 as build

##### METADATA #####
LABEL base.image="continuumio/miniconda3:4.7.10"
LABEL software="ImageMagick"
LABEL software.version="7.0.8"
LABEL software.description="Image manipulation suite"
LABEL software.website="https://imagemagick.org/"
LABEL software.documentation="https://imagemagick.org/"
LABEL software.license="custom"
LABEL software.tags="Image"
LABEL maintainer="alexander.kanitz@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://github.com/ImageMagick/ImageMagick/blob/master/LICENSE"

##### VARIABLES #####
ENV PACKAGES wget ca-certificates

##### INSTALLATION #####
RUN apt-get update -y
RUN apt-get install -y --no-install-recommends ${PACKAGES}
RUN wget https://imagemagick.org/download/binaries/magick

#### CLEAN IMAGE ####
FROM ubuntu:18.04
ENV PACKAGES libfontconfig1 libx11-dev libharfbuzz-dev libfribidi-bin
COPY --from=build magick /usr/bin/magick_raw
RUN echo '#!/bin/bash\n/usr/bin/magick_raw --appimage-extract-and-run "$@"' > /usr/bin/magick && \
    chmod +x /usr/bin/magick_raw && \
    chmod +x /usr/bin/magick && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    rm -rf /var/lib/apt/lists/*


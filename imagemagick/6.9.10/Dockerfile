##### BASE IMAGE #####
FROM debian:buster-slim as build

##### METADATA #####
LABEL base.image="debian:buster-slim" \
      software="ImageMagick" \
      software.version="6.9.10-23" \
      software.description="Image manipulation suite" \
      software.website="https://imagemagick.org/" \
      software.documentation="https://imagemagick.org/" \
      software.license="https://github.com/ImageMagick/ImageMagick/blob/master/LICENSE" \
      software.tags="image conversion" \
      maintainer="alexander.kanitz@unibas.ch" \
      maintainer.organisation="Biozentrum, University of Basel" \
      maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland" \
      maintainer.lab="Zavolan Lab"

###### INSTALLATION #####
RUN apt-get update -y
RUN apt-get install -y imagemagick

##### CLEAN IMAGE ####
FROM debian:buster-slim
COPY --from=build /usr/bin/convert /usr/bin/convert
ENV PACKAGES libmagickcore-6.q16-6 libmagickwand-6.q16-6 ghostscript
RUN apt-get update -y && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    rm -rf /var/lib/apt/lists/*


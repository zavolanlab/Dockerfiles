##### BASE IMAGE #####
FROM alpine:3.7

##### METADATA #####
LABEL base.image="alpine:3.7"
LABEL version="1"
LABEL software="Alpine"
LABEL software.version="3.7"
LABEL software.description="Alpine image with basics like bash and wget"
LABEL software.website="https://alpinelinux.org/"
LABEL software.documentation="https://wiki.alpinelinux.org/wiki/Main_Page"
LABEL software.license="https://alpinelinux.org/posts/Alpine-3.7.0-released.html"
LABEL software.tags="General"
LABEL maintainer="christina.herrmann@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### INSTALL #####
RUN apk update \
  && apk add --no-cache bash gawk sed grep bc coreutils perl wget gzip

##### WORKING DIRECTORY #####
WORKDIR /data/

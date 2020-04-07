FROM rocker/r-ver:3.6.3 as build
MAINTAINER Alex Kanitz "alexander.kanitz@alumni.ethz.ch"

## Install system dependencies
ENV PACKAGES build-essential git libcurl4-gnutls-dev libjpeg-dev libpng-dev libxml2-dev python3 python3-pip python3-setuptools python3-dev python-configparser zlib1g-dev
RUN apt-get update --yes && apt-get install -y --no-install-recommends ${PACKAGES} \
  build-essential \
  git \
  libcurl4-gnutls-dev \
  libjpeg-dev \
  libpng-dev \
  libxml2-dev \
  python3 \
  python3-pip \
  python3-setuptools \
  python3-dev \
  python-configparser \
  zlib1g-dev

## Install R packages from GitHub
RUN git clone https://github.com/l-modolo/TEtools && \
    cd /TEtools && \
    git checkout d557c54acd112b1ed12cc69811f19a68fb97126a

## Start with clean image
FROM rocker/r-ver:3.6.3
COPY --from=build /TEtools /TEtools
ENV PATH "/TEtools:${PATH}"

# Install system dependencies
ENV PACKAGES libcurl4-gnutls-dev libjpeg-dev libpng-dev libxml2-dev zlib1g-dev python3
WORKDIR /TEtools
RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    chmod +x fix_rosette.py PingPong.py TEcount.py TEdiff.R && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf .git

## Install R & Bioconductor packages
RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("gplots","ggplot2", "RColorBrewer"))' && \
  R -e 'install.packages("BiocManager")' && \
  R -e 'BiocManager::install()' && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("gplots","ggplot2", "RColorBrewer"))' && \
  R -e 'BiocManager::install(c("DESeq2"))'

CMD ["TEcount.py"]

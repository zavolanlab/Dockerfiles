##### BASE IMAGE #####
FROM ubuntu:18.04

##### METADATA #####
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="Whippet.jl"
LABEL software.version="0.11"
LABEL software.description="Lightweight and Fast; RNA-seq quantification at the event-level"
LABEL software.website="https://github.com/timbitz/Whippet.jl"
LABEL software.documentation="https://github.com/timbitz/Whippet.jl"
LABEL software.license="https://github.com/timbitz/Whippet.jl/blob/master/LICENSE"
LABEL software.tags="Genomics,Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"


RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes make g++ libgslcblas0 libglib2.0-dev libgtk2.0 wget unzip

# install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz \
  && tar xzvf julia-0.6.4-linux-x86_64.tar.gz \
  && ln -s /julia-9d11f62bcb/bin/julia /usr/bin/julia 

# install Whippet.jl
RUN julia -e 'Pkg.add("Whippet")' \
  && julia -e 'using Whippet'

WORKDIR ~/.julia/v0.6/Whippet/bin


#  && rm /usr/local/bin/julia-0.6.4-linux-x86_64.tar.gz \
#  && apt-get remove --purge --yes make g++  \
#  && apt-get autoremove --purge --yes

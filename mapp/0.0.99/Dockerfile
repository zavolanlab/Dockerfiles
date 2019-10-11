###############################################################################
#
#   Dockerfile for the container
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Biozentrum_University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 06-10-2019
#   LICENSE: GPL
#
###############################################################################

# BASE IMAGE
FROM continuumio/miniconda3:4.7.10

# METADATA
LABEL base.image="miniconda3:4.7.10"
LABEL version="0.0.99"
LABEL software="MAPP"
LABEL software.version="0.0.99"
LABEL software.description=\
"Dependencies for the pipeline to analyse motif activity on alternative \
splicing and polyadenylation."
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation=\
"Biozentrum, University of Basel, Swiss Institute of Bioinformatics"
LABEL maintainer.location=\
"Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"

# COPY THE YAML & INSTALL SOFTWARE WITH CONDA
COPY packages.yaml .
RUN conda env update --name base --file packages.yaml
RUN conda clean -all
##### BASE IMAGE #####
FROM continuumio/miniconda3:4.7.10

##### METADATA #####
LABEL base.image="continuumio/miniconda3:4.7.10"
LABEL version="1.0"
LABEL software=""
LABEL software.version=""
LABEL software.description="R libraries for data analysis and statistical modeling"
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license=""
LABEL software.tags="Data Analysis Statistical Modeling"
LABEL maintainer="maciej.bak@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license=""

##### VARIABLES #####
# Use variables for convenient updates/re-usability

# COPY THE YAML & INSTALL SOFTWARE WITH CONDA
COPY environment.yaml .
RUN conda env update --name base --file environment.yaml
RUN conda clean -all

# VARIABLES
ARG WORKDIR="/home/container/"
ARG USER="USER"
ARG GROUP="GROUP"
ENV PATH="${WORKDIR}:${PATH}"

# CREATE USER
RUN groupadd -r ${GROUP} && useradd --no-log-init -r -g ${GROUP} ${USER}

# SET ENVIRONMENT
WORKDIR ${WORKDIR}
RUN chown -R ${USER}:${GROUP} ${WORKDIR} && chmod 700 ${WORKDIR}
USER ${USER}

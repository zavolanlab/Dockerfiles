FROM perl

##### METADATA #####
LABEL base.image="perl"
LABEL version="1"
LABEL software="ucsc_pairwise_alignments_custom_scripts_perl"
LABEL software.version="0.1"
LABEL software.description="Custom scripts for ucsc alignments (perl) originally written from Dimos Gaidatzis"
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license="link to license(file) of original software"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="foivos.gypas@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 0.1


COPY lib /app/lib/
COPY *.pl /usr/bin/

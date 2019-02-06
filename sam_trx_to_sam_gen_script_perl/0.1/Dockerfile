##### BASE IMAGE #####
FROM perl:5.28

##### METADATA #####
LABEL base.image="perl:5.28"
LABEL version="1"
LABEL software="sam_trx_to_gen_script_perl"
LABEL software.version="0.1"
LABEL software.description="Perl script sam_trx_to_gen originally written by Alexander Kanitz"
LABEL software.website=""
LABEL software.documentation=""
LABEL software.license="link to license(file) of original software"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="paula.iborradetoledo@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Klingelbergstrasse 50/70, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Zavolan Lab"
LABEL maintainer.license="https://spdx.org/licenses/Apache-2.0"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION 0.1


COPY *.pl /usr/bin/

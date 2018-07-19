##### BASE IMAGE #####
FROM {{cookiecutter.base_image}}:{{cookiecutter.base_tag}}

##### METADATA #####
LABEL base.image="{{cookiecutter.base_image}}:{{cookiecutter.base_tag}}"
LABEL version="{{cookiecutter.dockerfile_version}}"
LABEL software="{{cookiecutter.orig_software}}"
LABEL software.version="{{cookiecutter.orig_software_version}}"
LABEL software.description="{{cookiecutter.orig_software_description}}"
LABEL software.website="{{cookiecutter.orig_software_url}}"
LABEL software.documentation="{{cookiecutter.orig_software_documentation}}"
LABEL software.license="{{cookiecutter.orig_software_license}}"
LABEL software.tags="{{cookiecutter.orig_software_tag}}"
LABEL maintainer="{{cookiecutter.maintainer}}"
LABEL maintainer.organisation="{{cookiecutter.organization}}"
LABEL maintainer.location="{{cookiecutter.location}}"
LABEL maintainer.lab="{{cookiecutter.lab}}"
LABEL maintainer.license="{{cookiecutter.license}}"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION {{cookiecutter.orig_software_version}}

##### INSTALL #####
RUN apt-get update -y \
  && apt-get install -y wget python build-essential libz-dev\
# Do your install stuff, merge with '&&' whenever possible, don't create unnecessary layers
  && wget https://.../v${SOFTWARE_VERSION}/software-${SOFTWARE_VERSION}.tar.gz \
  && tar -zxvf software-${SOFTWARE_VERSION}.tar.gz \
  && cd software \
  && make \
# Copy binaries to location in PATH
  && cp bin/* /usr/bin \
# Remove files no longer needed
  && cd .. \
  && rm -r /intermediate/folder software-${SOFTWARE_VERSION}.tar.gz \
# Clean up
  && apt-get purge -y wget \
  && apt-get autoremove -y && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \

##### SET PATH (if not in /usr/bin/) #####
ENV PATH /path/to/binary:$PATH

##### WORKING DIRECTORY #####
WORKDIR /data/

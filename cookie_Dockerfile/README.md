# cookie_Dockerfile
cookiecutter template for Dockerfiles

Simple template creating the following folder structure:  
```
  /Software_name/  
    /Software_version/  
      - Dockerfile
```
You will be prompted for stuff like:
- base image and tag
- software name, version and metadata
- maintainer information

For comprehensive fields see 'cookiecutter.json'

# About
This Dockerfile template follows standards from [BioContainers](http://biocontainers.pro/docs/developer-manual/biocotainers-dockerfile/). Please also read their [best practices](http://biocontainers.pro/docs/developer-manual/best-practices/) on how to create efficient docker images. You'll also find additional links to docker-how-to pages there.

*Note on deviations from BioContainers standards:*
- [MAINTAINER deprecated](https://docs.docker.com/engine/reference/builder/#maintainer-deprecated), using `LABEL maintainer="me@home.org"` instead.
- additionally included `maintainer.organization, -location, -lab and -license`.


# How to
1. [Install cookiecutter](https://cookiecutter.readthedocs.io/en/latest/installation.html)
2. cd to destination folder (where you want to create your Dockerfile-folder)
3. `cookiecutter https://github.com/ninsch3000/cookie_Dockerfile`
4. Fill in the values
5. Edit and extend the created Dockerfile, keeping the [best practices](http://biocontainers.pro/docs/developer-manual/best-practices/) in mind.

If you want to change the default values (especially for maintainer info), clone the repo and edit `cookiecutter.json`. You can then also call cookiecutter on your local copy of cookie_Dockerfile.

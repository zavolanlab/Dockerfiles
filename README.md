[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build_status](https://travis-ci.com/zavolanlab/Dockerfiles.svg?branch=master)](https://travis-ci.com/zavolanlab/Dockerfiles)

# Zavolan lab Dockerfile repository

Contains recipes for building Docker images for bioinformatics tools.
The main focus is on scripts/tools built within the [Zavolan lab], but
third-party tools are also included, if no well-maintained official image
repository exists (or has existed at the time of creation) for these tools.

The corresponding images can be found on the lab's [Docker Hub space].

## Contributing

Your contributions are highly appreciated!

To add a Dockerfile for a new tool or a different version of an existing tool,
please follow these simple rules:

1. If you are already a member of the `zavolanlab` organization here on GitHub,
   please clone the repository and create a _feature branch_ off of branch
   `master` (see [Git Flow] for details on why/how to do that). If you are an
   external contributor, please fork the repository instead.
2. Please create a new directory with the tool's name and/or version. It is
   important that you **stick to the directory structure/naming conventions:**
   `<root_dir>/<tool_name>/<version>/Dockerfile`. If you are
   unsure, have a look at some of the existing examples.
3. In order to reduce image size and possible security risks, try to only add
   software that is required to run the tool and follow [best practices] for
   writing Dockerfiles. Do not forget to include relevant metadata as well
   (have a look at the existing examples).
4. Document your code and update all relevant documentation, if necessary. Add
   your name and GitHub profile URL to the [list of
   contributors](contributors.md).
5. Issue a pull request (PR). Upon filing, images are automatically built for
   all new and modified Dockerfiles within the scope of the PR. PRs are not
   merged into the main branch (`master`) unless this check passes and at least
   one member of the Zavolan lab approves of the PR. Once merged, images are
   automatically pushed to the lab's [Docker Hub space].

For more information (particularly if you have are new to Git/GitHub) or other
ways of how you can contribute, please take a look at the [contributor
instructions](CONTRIBUTING.md).

## License

This project is covered by the [Apache License 2.0] also [shipped with this
repository](LICENSE).

[Apache License 2.0]: <https://www.apache.org/licenses/LICENSE-2.0>
[best practices]: <https://docs.docker.com/develop/develop-images/dockerfile_best-practices/>
[Docker Hub space]: <https://hub.docker.com/orgs/zavolab/repositories>
[Git Flow]: <https://datasift.github.io/gitflow/IntroducingGitFlow.html>
[Zavolan Lab]: <https://zavolan.biozentrum.unibas.ch/>


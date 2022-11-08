#!/bin/bash
# Build and export docker image as a *.tar.gz
# The image can then be loaded using "docker load --input swig_dev_MMDDYYYY.tar.gz"


# Default name to "swig_dev"
if [ -z "${VARIABLE}" ]; then
    DOCKER_NAME='swig_dev'
else
    DOCKER_NAME=${1}
fi

# Add today's date
DOCKER_NAME+="_$(date +'%m%d%Y')"

# Remove image if it already exists
docker rmi ${DOCKER_NAME}

# Remove tarball if it already exists
TAR_GZ_FILE="${DOCKER_NAME}.tar.gz"
[ -e ${TAR_GZ_FILE} ] && rm ${TAR_GZ_FILE}

# Build docker image and save image file as gzipped tarball (*.tar.gz)
docker buildx build --no-cache --tag=${DOCKER_NAME} . && \
docker save ${DOCKER_NAME}:latest | gzip > ${TAR_GZ_FILE}

#!/bin/bash

# Build and run docker image
docker buildx build --tag=swig_dev -f docker/Dockerfile . && \
docker run --name=swig_mrd_server --publish=9002:9002 --rm -it swig_dev

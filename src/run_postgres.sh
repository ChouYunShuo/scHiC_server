#!/bin/bash

docker run -it --name postgres \
    -e POSTGRES_DB=docker_hic \
    -e POSTGRES_USER=yunshuoc \
    -e POSTGRES_PASSWORD=ilovecoding123 \
    -e POSTGRES_HOST_AUTH_METHOD=trust \
    -p 5434:5432 \
    postgres:16-alpine
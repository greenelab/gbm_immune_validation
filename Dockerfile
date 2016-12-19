FROM ubuntu:16.04
MAINTAINER "Greg Way" <gregway@mail.med.upenn.edu>

RUN apt-get update && apt-get install -y --no-install-recommends \
            libcurl4-openssl-dev \
            libxml2-dev \
            ed \
            git \
            vim \
            curl \
            r-base=3.2.3-4 \
            r-base-dev=3.2.3-4 \
            r-recommended=3.2.3-4

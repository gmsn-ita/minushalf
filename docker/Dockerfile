FROM python:3.14-slim-bookworm

ARG DEBIAN_FRONTEND=noninteractive
ARG QE_VERSION=7.5

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    libopenblas-dev \
    liblapack-dev \
    libfftw3-dev \
    ca-certificates \
    curl \
    git \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN curl -L \
    "https://gitlab.com/QEF/q-e/-/archive/qe-${QE_VERSION}/q-e-qe-${QE_VERSION}.tar.gz" \
    -o qe.tar.gz \
    && tar -xzf qe.tar.gz \
    && rm qe.tar.gz \
    && mv "q-e-qe-${QE_VERSION}" quantum-espresso

WORKDIR /opt/quantum-espresso

RUN ./configure \
    && make -j2 all

ENV PATH="/opt/quantum-espresso/bin:/opt/quantum-espresso/upflib:${PATH}"

WORKDIR /workspace

CMD ["/bin/bash"]

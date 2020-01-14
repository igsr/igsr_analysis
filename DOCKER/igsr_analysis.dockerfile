FROM python:3.7-slim-buster

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build the image used in the different IGSR tasks"

# Install packages
RUN apt-get update \
 && apt-get -y --no-install-recommends install \
    build-essential \
    git \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
 && pip install --upgrade pip && \
    pip install --no-cache-dir pandas

# cloning igsr-analysis
WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis
ENV PATH=/bin/:/lib/igsr_analysis/scripts/VCF/QC/BENCHMARKING_TRUESET:${PATH}

#install HTSlib
WORKDIR /tmp/
RUN git clone --branch 1.9 https://github.com/samtools/htslib.git
WORKDIR /tmp/htslib/
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
#RUN rm -r /tmp/htslib


#install BCFTools
WORKDIR /tmp/
RUN git clone --branch 1.9 https://github.com/samtools/bcftools.git
WORKDIR /tmp/bcftools/
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
#RUN rm -r /tmp/bcftools

# getting rid of installed dependencies
RUN apt-get remove --purge -y git

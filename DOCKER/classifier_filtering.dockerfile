#parent image
FROM ubuntu:latest

MAINTAINER Ernesto Lowy <ernesto@ebi.ac.uk>

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install git \
    				   build-essential \
				   autoconf \
				   zlib1g-dev \
				   libbz2-dev \
				   liblzma-dev \
				   libhts-dev  \
				   python3 \
				   python3-pip \
				   && apt-get clean
				   
WORKDIR tmp/

#prepare Python
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/pip3 /usr/bin/pip

#install HTSlib
RUN git clone https://github.com/samtools/htslib.git
WORKDIR htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
WORKDIR tmp/bcftools

#install BCFTools
RUN git clone https://github.com/samtools/bcftools.git
WORKDIR bcftools/
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install

#install igsr-analysis libraries
WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis

RUN pip install pandas sklearn
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
				   libcurl4-openssl-dev \
				   libssl-dev \
				   && apt-get clean
				   
WORKDIR tmp/

#prepare Python
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/pip3 /usr/bin/pip

#install HTSlib
RUN git clone --branch 1.9 https://github.com/samtools/htslib.git
WORKDIR htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
WORKDIR /tmp/bcftools
RUN rm -r /tmp/htslib

#install BCFTools
RUN git clone --branch 1.9 https://github.com/samtools/bcftools.git
WORKDIR bcftools/
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
WORKDIR /tmp/vcflib
RUN rm -r /tmp/bcftools

#install vcflib
RUN git clone --recursive https://github.com/vcflib/vcflib.git
WORKDIR vcflib/
RUN make
RUN cp bin/vcfallelicprimitives /bin/
 
#install igsr-analysis libraries
WORKDIR /lib
RUN rm -r /tmp/vcflib/
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis
ENV PATH=/bin/:${PATH}

#install vt
WORKDIR /tmp/vt
RUN git clone https://github.com/atks/vt.git
WORKDIR vt/
RUN make
RUN cp vt /bin/
WORKDIR /root/
RUN rm -r /tmp/vt

RUN pip install pandas sklearn
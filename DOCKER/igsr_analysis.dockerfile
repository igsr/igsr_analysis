#parent image
FROM alpine:3.7

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build the image used in the different IGSR tasks" 

RUN apk add --no-cache git

# cloning igsr-analysis
WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis
ENV PATH=/bin/:/lib/igsr_analysis/scripts/VCF/QC/BENCHMARKING_TRUESET:${PATH}
				   

FROM alpine:3.14 AS builder
LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build run igsr-analysis unit tests" 

RUN apk add --no-cache git build-base zlib-dev bzip2-dev xz-dev curl-dev wget
WORKDIR tmp/

#install VG
WORKDIR vg
RUN wget https://github.com/vgteam/vg/releases/download/v1.34.0/vg && chmod +x vg

FROM frolvlad/alpine-miniconda3
RUN apk --no-cache add ca-certificates git zlib-dev curl-dev xz-dev bzip2-dev
WORKDIR /usr/local/bin/
COPY --from=builder /tmp/vg ./
RUN conda install -c bioconda gatk bcftools bedtools

WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis/
# Run tests
RUN pip install pytest
WORKDIR /lib/igsr_analysis/tests
RUN ["pytest", "-v", "--junitxml=reports/result.xml", "test_BCFTools.py","test_VgToolkit.py", "test_GATK_HC.py", "test_GATK_HC.py", "test_BEDTools.py"]

CMD tail -f /dev/null

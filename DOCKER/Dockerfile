FROM alpine:3.14 AS builder
LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build run igsr-analysis unit tests" 

RUN apk add --no-cache git build-base zlib-dev bzip2-dev xz-dev curl-dev wget
WORKDIR tmp/
# install BCFTools
RUN git clone --recurse-submodules git://github.com/samtools/htslib.git && git clone git://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
WORKDIR /tmp/
RUN rm -rf htslib && rm -rf bcftools

#install VG
WORKDIR vg
RUN wget https://github.com/vgteam/vg/releases/download/v1.34.0/vg && chmod +x vg

FROM alpine:3.14  
RUN apk --no-cache add ca-certificates git zlib-dev curl-dev xz-dev bzip2-dev py3-pip
RUN pip install --no-cache-dir pytest
WORKDIR /usr/local/bin/
COPY --from=builder /usr/local/bin/bcftools ./
COPY --from=builder /tmp/vg ./
WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis/
# Run tests
WORKDIR /lib/igsr_analysis/tests
RUN ["pytest", "-v", "--junitxml=reports/result.xml", "test_BCFTools.py"]

CMD tail -f /dev/null

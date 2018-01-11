FROM ubuntu:16.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "Camoco is a fully-fledged software package for building co-expression networks and analyzing the overlap interactions among genes."

RUN apt-get -y update && apt-get install -y \
        curl \
        wget \
        git \
        gcc \
        build-essential \
        libqt5gui5 \
        python3

RUN mkdir -p /src/ && \
    mkdir -p /data/ && \
    mkdir -p /cobdata/ && \
    cd /src/ && \
    git clone https://github.com/LinkageIO/Camoco.git && \
    cd Camoco/ && \
    ./install.sh

WORKDIR /root

#ADD /home/rob/.camoco/databases /root/.camoco/databases

ENV LD_LIBRARY_PATH=/root/.camoco/lib:$LD_LIBRARY_PATH \ 
    PATH=/root/.camoco/bin:/root/.camoco/conda/envs/camoco/bin:$PATH


RUN [ "/bin/bash", "-c", "source activate camoco" ]

ENTRYPOINT [ "/root/.camoco/conda/envs/camoco/bin/camoco" ]
#ENTRYPOINT [ "/bin/bash" ]

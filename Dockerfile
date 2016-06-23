FROM ubuntu:15.10

RUN apt-get -y update && apt-get install -y \
        curl \
        wget \
        git \
        vim \
        gcc \
        build-essential \
        libqt5gui5 \
        python3

RUN mkdir -p /src/ && \
    mkdir -p /data/ && \
    mkdir -p /cobdata/ && \
    cd /src/ && \
    git clone https://github.com/schae234/Camoco.git && \
    cd Camoco/ && \
    ./install.sh

ENV LD_LIBRARY_PATH=/root/.camoco/lib:$LD_LIBRARY_PATH \ 
    PATH=/root/.camoco/bin:/root/.camoco/conda/bin:$PATH

ENTRYPOINT "/bin/bash"

FROM ubuntu:16.04
MAINTAINER Rob<rob@linkage.io>
LABEL Description "Camoco is a fully-fledged software package for building co-expression networks and analyzing the overlap interactions among genes."

RUN apt-get -y update && apt-get install -y \
    curl \
	lsb-release \ 
    wget \
    git \
    gcc \
    build-essential \
    libqt5gui5 \
	apt-transport-https \
    python3

RUN mkdir -p /src/ && \
    mkdir -p /data/ && \
    mkdir -p /cobdata/ && \
    cd /src/ && \
    git clone https://github.com/LinkageIO/Camoco.git && \
    cd Camoco/ && \
    ./install.sh

# Activate the camoco env
ENV PATH=$HOME/.camoco/bin/:$PATH
RUN /bin/bash -c "source $HOME/.camoco/conda/bin/activate camoco"

ENV LD_LIBRARY_PATH=/root/.camoco/lib:$LD_LIBRARY_PATH 
ENV PATH=/root/.camoco/bin:/root/.camoco/conda/envs/camoco/bin:$PATH

RUN cd /src/Camoco

# Build a refgen
RUN camoco build-refgen \
    /src/Camoco/tests/raw/RefGen/ZmB73_5b_FGS.gff.gz \
    "Zm5bFGS" "Maize FGS" 5b "Zea mays"
# Build a COB
RUN camoco build-cob \ 
    /src/Camoco/tests/raw/Expr/RNASEQ/MaizeRNASeqTissue.tsv.bz2 \
    ZmTissue "Sekhon et el Tissue Atlas" \
    Zm5bFGS --rawtype RNASEQ

#ENTRYPOINT [ "/bin/bash" ]
ENTRYPOINT [ "/root/.camoco/conda/envs/camoco/bin/camoco" ]

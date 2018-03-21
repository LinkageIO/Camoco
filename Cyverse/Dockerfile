FROM ubuntu:14.04
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

# Install iRods
RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add -
RUN echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" | tee /etc/apt/sources.list.d/renci-irods.list
RUN apt-get update
RUN apt-get -y install irods-server irods-database-plugin-postgres

ENV IRODS_HOST=data.iplantcollaborative.org
ENV IRODS_PORT=1247
ENV IRODS_USER_NAME=anonymous
ENV IRODS_ZONE_NAME=iplant

RUN mkdir -p /root/.irods/
RUN echo '{"irods_zone_name": "iplant","irods_user_name": "anonymous","irods_host": "data.iplantcollaborative.org","irods_port": 1247}' | tee /root/.irods/irods_environment.json

#RUN apt-get update ; apt-get upgrade -y ; apt-get install -y wget libcurl4-gnutls-dev
#RUN wget ftp://ftp.renci.org/pub/irods/releases/4.1.3/ubuntu14/irods-icommands-4.1.3-ubuntu14-x86_64.deb -O /tmp/icommands.deb
#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y `dpkg -I /tmp/icommands.deb | sed -n 's/^ Depends: //p' | sed 's/,//g'`
#RUN dpkg -i /tmp/icommands.deb
     
#RUN cd /src && \
#    git clone https://github.com/irods/irods_client_fuse && \
    

WORKDIR /root

#ADD /home/rob/.camoco/databases /root/.camoco/databases

ENV LD_LIBRARY_PATH=/root/.camoco/lib:$LD_LIBRARY_PATH \ 
    PATH=/root/.camoco/bin:/root/.camoco/conda/envs/camoco/bin:$PATH

#RUN [ "/bin/bash", "-c", "source activate camoco" ]

ENTRYPOINT [ "/bin/bash" ]
#ENTRYPOINT [ "/root/.camoco/conda/envs/camoco/bin/camoco" ]

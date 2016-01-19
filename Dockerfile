FROM ubuntu:15.10

RUN apt-get -y update && apt-get install -y \
        curl \
        wget \
        mcl \
        git \
	vim \
        python3 \
        python3-pip \
        python3-setuptools \
        cython3 \
        python3-nose \
        python3-numpy \
        python3-scipy \
        python3-pandas \
        python3-matplotlib \
        python3-numexpr \
        python3-apsw \
        python3-tables \
        python3-yaml \
        ipython3 \
        python3-patsy \
        python3-flask \
        python3-networkx \
        python3-termcolor \
        python3-pytest \
        python3-ipdb

RUN pip3 install statsmodels
RUN pip3 install --no-deps -U pandas

RUN echo '--- # YAML Camoco Configuration File \n\
logging: \n\
  log_level: verbose \n\
options: \n\
  alpha: 0.0001 \n\
  basedir: /cobdata/ \n\
  debug: true \n\
  testdir: /src/Camoco/tests/ \n\
test: \n\
  cob: NewRoot \n\
  force: \n\
    COB: true \n\
    Ontology: true \n\
    RefGen: true \n\
  gene: GRMZM2G000014 \n\
  num: 50 \n\
  ontology: ZmIonome \n\
  refgen: Zm5bFGS \n\
  term: Fe57' > ~/.camoco.conf && \
	mkdir -p /src/ && \
	mkdir -p /data/ && \
	mkdir -p /cobdata/ && \
        cd /src/ && \
        git clone https://github.com/monprin/Camoco.git && \
        cd Camoco && \
        python3 setup.py install && \
        python3 -c 'import camoco as co;'

VOLUME /cobdata

RUN git config --global credential.helper cache

RUN cd /home/ && \
    wget https://github.com/github/git-lfs/releases/download/v0.5.4/git-lfs-linux-amd64-0.5.4.tar.gz && \
    tar xzf git-lfs-linux-amd64-0.5.4.tar.gz && \
    rm -rf git-lfs-linux-amd64-0.5.4.tar.gz && \
    cd git-lfs-0.5.4/ && \
    mv git-lfs /bin/ && \
    cd /home/ && \
    rm -rf git-lfs-0.5.4/ && \
    git lfs init

ENTRYPOINT ["/bin/bash"]

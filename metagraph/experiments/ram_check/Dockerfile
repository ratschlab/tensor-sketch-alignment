FROM --platform=linux/x86_64 ubuntu:20.04 AS metagraph_dev_env

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    binutils-dev \
    ccache \
    curl \
    g++-10 \
    git \
    libboost-all-dev \
    libbrotli-dev \
    libbz2-dev \
    libdouble-conversion-dev \
    libevent-dev \
    libgflags-dev \
    libgoogle-glog-dev \
    libiberty-dev \
    libjemalloc-dev \
    liblz4-dev \
    liblzma-dev \
    libsnappy-dev \
    libssl-dev \
    libtool \
    libunwind-dev \
    make \
    pkg-config \
    python3 \
    python3-pip \
    python3-venv \
    vim \
    wget \
    zlib1g-dev

RUN pip3 install parameterized==0.7.1 "cmake>=3.19"

ENV LD_LIBRARY_PATH=/usr/local/lib
ENV CC /usr/bin/gcc-10
ENV CXX /usr/bin/g++-10


RUN apt-get update && apt-get install -y \
    libatomic1 \
    libcurl4-nss-dev \
    libgomp1 \
    libhts-dev \
    libjemalloc2 \
    python3 \
    python3-pip 

RUN git clone https://github.com/alexandrumeterez/metagraph.git &&\
    cd metagraph &&\
    git checkout speed-increase &&\
    git pull &&\
    git submodule update --init --recursive &&\
    git submodule sync &&\
    git submodule update --init --recursive &&\
    oldpath=$PWD &&\
    cd metagraph/external-libraries/sdsl-lite &&\
    ./install.sh $PWD &&\
    cd $oldpath &&\
    mkdir metagraph/build &&\
    cd metagraph/build &&\
    cmake .. &&\
    make -j 7 


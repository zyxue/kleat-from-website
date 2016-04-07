# https://raw.githubusercontent.com/sjackman/docker-linuxbrew/master/linuxbrew-core/Dockerfile

FROM ubuntu
MAINTAINER Shaun Jackman <sjackman@gmail.com>

RUN localedef -i en_US -f UTF-8 en_US.UTF-8

RUN apt-get update
RUN apt-get install -y curl g++ gawk m4 make patch ruby tcl

RUN useradd -m -s /bin/bash linuxbrew
RUN echo 'linuxbrew ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER linuxbrew
WORKDIR /home/linuxbrew
ENV PATH /home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:$PATH
ENV SHELL /bin/bash
RUN yes |ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/linuxbrew/go/install)"
RUN brew doctor || true

#####################above are about setting up linuxbrew######################

# install transabyss requirements
RUN sudo apt-get install -y python-pip python-dev zlib1g-dev libxml2-dev

RUN sudo pip install pysam
RUN sudo pip install biopython
RUN sudo pip install python-igraph

RUN sudo apt-get install -y git
RUN sudo apt-get install -y default-jre
RUN sudo apt-get install -y libncurses5 libncurses5-dev

RUN sudo apt-get -yf install \
automake \
bsdmainutils \
build-essential \
curl \
gcc \
libboost-all-dev \
libopenmpi-dev \
libpython-dev \
libsparsehash-dev \
libsqlite3-dev \
openmpi-bin \
python-dev \
python-setuptools \
ssh \
unzip \
zlib1g-dev

RUN sudo apt-get install -y libopenblas-dev libopenblas-base liblapack3 liblapack-dev liblapacke liblapacke-dev

# by default installing bowtie2 fails without --without-tbb
RUN brew install -dv homebrew/science/bowtie2 --without-tbb

RUN brew install -dv homebrew/science/trans-abyss

# trans-abyss doesn't work with abyss-1.9.0 selected by linuxbrew, so install
# an old version instead
RUN brew uninstall abyss

RUN curl -OL "https://github.com/bcgsc/abyss/releases/download/1.5.2/abyss-1.5.2.tar.gz" \
    && tar zxf abyss-1.5.2.tar.gz \
    && cd abyss-1.5.2 \
    && ./configure --with-mpi=/usr/lib/openmpi \
    && make -j 8 \
    && sudo make install \
    && cd .. && rm -rfv abyss-1.5.2*

RUN brew install -dv homebrew/science/biobloomtools

RUN brew install -dv homebrew/science/bwa

RUN sudo apt-get install -y python-pip \
    && sudo pip install -U colorlog ruffus

ENV HOME /home/linuxbrew

# use an older version of gmap-gsnap
ENV PATH=${HOME}/local/gmap/bin:${PATH}
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-03-28.v2.tar.gz  \
    && tar zxf gmap-gsnap-2014-03-28.v2.tar.gz \
    && cd gmap-2014-03-28 \
    && ./configure --prefix=${HOME}/local/gmap \
    && make \
    && make install \
    && cd .. && rm -rf gmap-gsnap-2014-03-28.v2.tar.gz gmap-2014-03-28

ENV PROJECT_ID=tasrkleat
ENV _PROJECT_DIR ${HOME}/local/${PROJECT_ID}
ENV PATH=${_PROJECT_DIR}:${PATH}
RUN mkdir -p ${_PROJECT_DIR}
ADD app/*.py ${_PROJECT_DIR}/

FROM debian:jessie

# TODO: install only the actually required boost libraries

RUN apt-get update -y && DEBIAN_FRONTEND=noninteractive apt-get install -y \
  libboost-all-dev \
  make \
  g++ \
  libreadline-dev \
  zlib1g-dev \
  libncurses5-dev


ADD . /mothur/

# for Linux TARGET_ARCH has to be commented
RUN sed -i -e 's/TARGET_ARCH /#TARGET_ARCH /' \
  -e 's/zlib.a/libz.a/' \
  -e 's/BOOST_LIBRARY_DIR=\".*\"/BOOST_LIBRARY_DIR=\"\/usr\/lib\/x86_64-linux-gnu\/\"/' \
  /mothur/makefile

RUN cd /mothur && make

ENV PATH "/mothur/:$PATH"

# example
# git clone -b v1.36.1 https://github.com/mothur/mothur.git
# cd mothur
# docker build -t mothur:v1.36.1 .


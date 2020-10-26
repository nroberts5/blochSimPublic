FROM ubuntu:xenial
LABEL Name=ubuntu_docker Version=0.0.1

# Update apps on the base image
RUN apt -y update && apt install -y \
    build-essential \
    gawk \
    xutils-dev \
    csh \
    tcsh \
    wget \
    git
RUN ln -sf /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/libz.so && ln -s /usr/bin/awk /bin/awk

# Install Newer Cmake
WORKDIR /usr/src/cmake
RUN apt -y remove --purge --auto-remove cmake && version=3.16 && build=2 && wget https://cmake.org/files/v$version/cmake-$version.$build-Linux-x86_64.sh --no-check-certificate
RUN mkdir /opt/cmake && version=3.16 && build=2 && printf 'y\nn\n' | sh cmake-$version.$build-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

RUN apt-get -y update && apt-get -y install \
    liblapack-dev \
    libblas-dev \
    libboost-dev \
    libarmadillo-dev
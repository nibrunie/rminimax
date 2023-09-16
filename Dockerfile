
# Pull base image (use Wily for now).
FROM ubuntu:22.04

# Set the maintainer
MAINTAINER Nicolas Brunie

# Install some base tools that we will need to get the risc-v
# toolchain working.
RUN apt update && apt install -y git \
        build-essential \
        autoconf \
	automake \
	autotools-dev \
	curl \
	libmpc-dev \
	libmpfr-dev \
	libgmp-dev \
	libusb-1.0-0-dev \
	gawk \
	build-essential \
	bison \
	flex \
	texinfo \
	gperf \
	libtool \
	patchutils \
	bc \
	zlib1g-dev \
	device-tree-compiler \
	pkg-config \
	libexpat-dev


# project specific dependencies
RUN apt update && apt install -y libfplll-dev gnuplot cmake
RUN apt update && apt install -y libflint-arb-dev libeigen3-dev 
#RUN apt install -y  arb libeigen-dev flint

# Make a working folder and set the necessary environment variables.
ENV LOCAL_INSTALL /opt/local
ENV NUMJOBS 1
RUN mkdir -p $LOCAL_INSTALL
RUN mkdir -p $LOCAL_INSTALL/include/

# Make a directory to download sources and build dependencies
ENV DEP_SRCS /home/dep_srcs/
RUN mkdir -p $DEP_SRCS

ENV PATH $LOCAL_INSTALL/bin:$PATH

WORKDIR $DEP_SRCS

# qsopt_ex
RUN git clone https://github.com/jonls/qsopt-ex.git
WORKDIR $DEP_SRCS/qsopt-ex
RUN ./bootstrap
RUN mkdir build
WORKDIR $DEP_SRCS/qsopt-ex/build
RUN ../configure --prefix=$LOCAL_SINTALL
RUN make && make install


# mpreal
WORKDIR $DEP_SRCS
RUN git clone https://github.com/advanpix/mpreal.git
RUN ln -s $DEP_SRCS/mpreal/mpreal.h /opt/local/include/

# building rminimax (code is binded to host directory)
WORKDIR /home/
RUN git clone https://github.com/nibrunie/rminimax.git
WORKDIR /home/rminimax
RUN mkdir build
WORKDIR /home/rminimax/build/
RUN cmake -DMPREAL_INCLUDE_DIR=/opt/local/include/ ..
RUN make

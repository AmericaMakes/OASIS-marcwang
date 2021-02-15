from continuumio/miniconda3:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gdb

COPY ./thirdparty/catch2 /tmp/catch2
RUN cd /tmp/catch2 
RUN cd /tmp/catch2 && \
    cmake -Bbin -H. -DBUILD_TESTING=OFF && \
    cd bin &&\
    make . install -j4

RUN rm -rf /tmp/catch2

COPY ./thirdparty/geogram /tmp/geogram
RUN cd /tmp/geogram && \
    sh -f configure.sh && \
    cd build/Linux64-gcc-dynamic-Release && \
    make . install -j4

RUN rm -rf /tmp/geogram

COPY ./dependencies.yml tmp/
RUN conda env update --file "tmp/dependencies.yml"
RUN rm -f tmp/dependencies.yml
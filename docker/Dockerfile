FROM python:3.8-slim AS build
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --no-install-recommends --no-install-suggests -yq \
    g++ make cmake git ca-certificates \
    libhdf5-serial-dev h5utils cmake cmake-curses-gui libboost-all-dev \
    libboost-all-dev xsdcxx libxerces-c-dev libhdf5-serial-dev h5utils \
    hdf5-tools libtinyxml-dev libxml2-dev libxslt1-dev libfftw3-dev \
    liblapacke-dev libpng-dev libopenblas-dev gfortran ninja-build && \
    rm -rf /var/lib/apt/lists/*

RUN  mkdir -p /opt/code

# ISMRMRD library
RUN cd /opt/code && \
    git clone --no-checkout https://github.com/ismrmrd/ismrmrd && \
    cd ismrmrd && \
    git checkout v1.8.0 && \
    mkdir build && \
    cd build && \
    cmake -G Ninja ../ && \
    ninja && \
    ninja install && \
    cd /usr/local/lib && \
    tar -czvf libismrmrd.tar.gz libismrmrd*

# BART
RUN cd /opt/code && \
    git clone --no-checkout https://github.com/mrirecon/bart && \
    cd bart && \
    git checkout v0.8.00 && \
    make SLINK=1 OPENBLAS=1 OMP=1 FFTW_WISDOM=1 -j $(nproc)

# Clone from GitHub
RUN cd /opt/code && \
    git clone https://github.com/ismrmrd/ismrmrd-python-tools.git && \
    git clone https://github.com/kspaceKelvin/python-ismrmrd-server.git && \
    python3 -m pip wheel --wheel-dir /wheels /opt/code/ismrmrd-python-tools/ && \
    rm -rf /opt/code/ismrmrd-python-tools

# Build remaining wheels
RUN python3 -m pip wheel --wheel-dir /wheels pyxb h5py numpy scipy ismrmrd

# ----- Start a clean build without build dependencies -----
FROM python:3.8-slim
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    TOOLBOX_PATH=/opt/code/bart

# Copy from last stage...
COPY --from=build /usr/local/lib/libismrmrd.tar.gz /usr/local/lib/
COPY --from=build /opt/code /opt/code
COPY --from=build /wheels /wheels

# Add libismrmrd
RUN cd /usr/local/lib && \
    tar -zxvf libismrmrd.tar.gz && \
    rm -f libismrmrd.tar.gz && \
    ldconfig

RUN python3 -m pip install -U pip && \
    python3 -m pip --no-cache-dir install --no-index -f /wheels pyxb h5py numpy scipy ismrmrd ismrmrd-python-tools

# Add local code
COPY . /opt/code/python-ismrmrd-server

CMD ["python3", "/opt/code/python-ismrmrd-server/main.py", "-H=0.0.0.0", "-p=9002", "-l=/tmp/fire.log"]

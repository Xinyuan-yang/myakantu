FROM registry.gitlab.com/epfl-lsms/cicd-templates/cicd-templates-cpp:trixie-slim

RUN apt -qq update && apt -qq -y install \
    g++ gfortran clang cmake ccache git \
    openmpi-bin libmumps-dev libscotch-dev \
    libboost-dev libopenblas-dev \
    libeigen3-dev libhdf5-dev libhdf5-openmpi-dev \
    python3 python3-dev python3-numpy python3-scipy python3-mpi4py \
    gmsh \
    && rm -rf /var/lib/apt/lists/*

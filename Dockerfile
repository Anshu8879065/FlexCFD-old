# Multi-stage Dockerfile for building and packaging the FlexCFD project
# - Builder stage uses Ubuntu 22.04, installs Homebrew (Linux) and PETSc via brew
# - Builder performs the full build and installs artifacts into /install
# - Final runtime stage is slim and only copies runtime artifacts (PETSc, app)

# --------------------
# Builder stage
# --------------------
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/home/builder \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

ARG USER=builder
ARG UID=1000

# Install only build-essential/tools and runtime libraries required for building
# NOTE: we intentionally do NOT install system MPI/OpenMPI or BLAS/LAPACK here;
# Homebrew will provide PETSc and its dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake ninja-build git curl ca-certificates \
    python3 python3-pip python3-venv pkg-config wget sudo \
    locales tzdata libssl-dev && \
    locale-gen en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

# Create a non-root user to run Homebrew and build the project
RUN useradd -m -d /home/${USER} -s /bin/bash -u ${UID} ${USER} && \
    echo "${USER} ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${USER} && \
    chmod 0440 /etc/sudoers.d/${USER}

USER ${USER}
WORKDIR /home/${USER}

# Ensure Homebrew paths are available in subsequent RUN layers
ENV PATH=/home/${USER}/.local/bin:/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:${PATH}

# Install Homebrew (non-interactive) for Linux and set up environment
RUN export NONINTERACTIVE=1 && \
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" && \
    echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' >> /home/${USER}/.profile && \
    eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"

# Use Homebrew to install PETSc (let brew bring in MPI and other dependencies)
RUN brew update && brew install petsc && brew install gcc@14 || true

# Make sure wrappers that expect e.g. gcc-14 can find a compiler (OpenMPI wrappers may hardcode a compiler name)
USER root
RUN ln -sf /home/linuxbrew/.linuxbrew/bin/gcc-14 /usr/bin/gcc-14 || true && \
    ln -sf /home/linuxbrew/.linuxbrew/bin/g++-14 /usr/bin/g++-14 || true && \
    ln -sf /home/linuxbrew/.linuxbrew/bin/gfortran-14 /usr/bin/gfortran-14 || true
USER ${USER}

# Set PETSc environment for the build stage
ENV PETSC_DIR=/home/linuxbrew/.linuxbrew/opt/petsc
ENV PETSC_ARCH=
ENV LD_LIBRARY_PATH=${PETSC_DIR}/lib:${LD_LIBRARY_PATH}

# Install Python tooling (conan)
RUN python3 -m pip install --user conan

# Copy the project into the builder
USER root
WORKDIR /home/${USER}
COPY --chown=${USER}:${USER} . /home/${USER}/app
USER ${USER}
WORKDIR /home/${USER}/app

# Build the project using conan + cmake + ninja and install into /install
RUN export PATH=/home/${USER}/.local/bin:$PATH && \
    export CC=/home/linuxbrew/.linuxbrew/bin/gcc-14 CXX=/home/linuxbrew/.linuxbrew/bin/g++-14 FC=/home/linuxbrew/.linuxbrew/bin/gfortran-14 && \
    export OMPI_CC=/home/linuxbrew/.linuxbrew/bin/gcc-14 OMPI_CXX=/home/linuxbrew/.linuxbrew/bin/g++-14 && \
    rm -rf build && \
    mkdir -p build && \
    cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/${USER}/install \
        -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_Fortran_COMPILER=${FC} && \
    cmake --build build -j"$(nproc)" && \
    cmake --install build --prefix /home/${USER}/install

# --------------------
# Runtime stage
# --------------------
FROM ubuntu:22.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Install minimal runtime packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# Create non-root runtime user
ARG USER=builder
ARG UID=1000
RUN useradd -m -d /home/${USER} -s /bin/bash -u ${UID} ${USER} || true

# Copy installed application and PETSc runtime artifacts from builder
# We copy the /install tree (cmake install prefix) and the Homebrew prefix so PETSc and MPI libs are available
COPY --from=builder /home/builder/install /usr/local
COPY --from=builder /home/linuxbrew/.linuxbrew /home/linuxbrew/.linuxbrew

# Fix permissions
RUN chown -R ${USER}:${USER} /home/${USER} /usr/local || true

# Set environment to find PETSc and brew-provided MPI
ENV PATH=/usr/local/bin:/home/linuxbrew/.linuxbrew/bin:/home/linuxbrew/.linuxbrew/sbin:${PATH}
ENV PETSC_DIR=/home/linuxbrew/.linuxbrew/opt/petsc
ENV LD_LIBRARY_PATH=/home/linuxbrew/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH}

USER ${USER}
WORKDIR /home/${USER}

# Default to interactive shell (override when running container)
CMD ["/bin/bash"]

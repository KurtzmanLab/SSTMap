FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

# Install system packages, including gcc-8 and g++-8
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    ca-certificates \
    curl \
    git \
    libgsl-dev \
    gcc-7 g++-7 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.2-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# Set working directory
WORKDIR /app/sstmap

RUN mkdir -p /app/output

# Copy SSTMap source code
COPY . .

# Use bash so conda works correctly
SHELL ["/bin/bash", "-lc"]

# Create conda environment
RUN conda env create -f environment.yml

# Make the conda environment active by default (NO conda activate needed)
ENV CONDA_PREFIX=/opt/conda/envs/sstmap
ENV PATH="$CONDA_PREFIX/bin:$PATH"

# Install any extra packages inside the conda environment
RUN pip install pyparsing==3.0.4

# Install SSTMap inside the conda environment
RUN python setup.py install


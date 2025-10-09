# Use a minimal base image
FROM ubuntu:22.04

# Set noninteractive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    git \
    build-essential \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Accept Anaconda Terms of Service
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Update conda and install mamba
RUN conda update -n base -c defaults conda -y && \
    conda install -n base -c conda-forge mamba -y

# Create bioinformatics environment with Python, R, Bioconductor, and tools
RUN mamba create -y -n bioinfo \
    -c conda-forge -c bioconda -c defaults \
    star=2.7.10a \
    samtools=1.13 \
    bedtools=2.30.0 \
    gffcompare=0.11.2 \
    minimap2=2.22 \
    pbccs=6.0.0 \
    lima=2.0.0 \
    isoseq3=3.4.0 \
    deeptools=3.5.1 \
    pybigwig=0.3.18 \
    bx-python=0.8.11 \
    numpy=1.21.2 \
    pandas=1.3.2 \
    scipy=1.7.1 \
    scikit-learn=0.24.2 \
    joblib=1.0.1 \
    threadpoolctl=2.2.0 \
    r-base=4.3.1 \
    r-devtools \
    r-tidyverse \
    r-ggplot2 \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-genomicfeatures \
    bioconductor-genomicranges \
    bioconductor-rtracklayer \
    bioconductor-tximport \
    && conda clean -afy

# Install Python cupcake library via pip inside the bioinfo environment
RUN /opt/conda/envs/bioinfo/bin/pip install cDNA_Cupcake

# Activate environment by default
SHELL ["conda", "run", "-n", "bioinfo", "/bin/bash", "-c"]

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]
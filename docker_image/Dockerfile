# Use a minimal Ubuntu base
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------------------------
# 1. System dependencies
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    git \
    build-essential \
    curl \
    unzip \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# 2. Python virtual environment
# ------------------------------------------------------------------------------
RUN python3 -m venv /opt/bioinfo
ENV PATH="/opt/bioinfo/bin:$PATH"

# Upgrade pip and setuptools
RUN pip install --upgrade pip setuptools wheel

# ------------------------------------------------------------------------------
# 3. Python packages
# ------------------------------------------------------------------------------

RUN pip install \
    numpy==1.23.5 \
    pandas==1.5.3 \
    scipy==1.10.1 \
    scikit-learn==1.2.2 \
    joblib==1.3.2 \
    threadpoolctl==3.1.0 \
    pybigwig==0.3.18 \
    git+https://github.com/bxlab/bx-python.git \
    deeptools==3.5.1

# ------------------------------------------------------------------------------
# 4. Command-line bioinformatics tools
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz -O /tmp/star.tar.gz && \
    tar -xzf /tmp/star.tar.gz -C /tmp && \
    cd /tmp/STAR-2.7.10a/source && make STAR && \
    mv STAR /usr/local/bin/ && rm -rf /tmp/STAR-2.7.10a /tmp/star.tar.gz

# minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2 -O /tmp/minimap2.tar.bz2 && \
    tar -xjf /tmp/minimap2.tar.bz2 -C /usr/local/bin/ && \
    rm -rf /tmp/minimap2.tar.bz2

# ------------------------------------------------------------------------------
# 5. Working directory and entrypoint
# ------------------------------------------------------------------------------
WORKDIR /data
CMD ["bash"]
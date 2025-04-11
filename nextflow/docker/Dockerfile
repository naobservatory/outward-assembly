# Use a lightweight base (you can also pick a specific tag, e.g. ubuntu:22.04)
FROM ubuntu:latest

# Install zstd (for zstdcat) and dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        zstd \
        wget \
        openjdk-11-jre-headless && \
    rm -rf /var/lib/apt/lists/*

# Download and install a specific version of BBTools (BBMap suite)
# This is the latest version as of 2025-02-17
RUN mkdir -p /opt/bbtools && \
    cd /opt/bbtools && \
    wget -q https://sourceforge.net/projects/bbmap/files/BBMap_39.17.tar.gz/download -O bbmap.tar.gz && \
    tar -xzvf bbmap.tar.gz && \
    rm bbmap.tar.gz

# Add BBTools to the PATH
ENV PATH="/opt/bbtools/bbmap:$PATH"

# Set the working directory to /mnt
WORKDIR /mnt

# Combined Dockerfile with Python, pandas, gffcompare, gffutils, and gawk
FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive

# Update package lists and install system dependencies
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    wget \
    ca-certificates \
    libgomp1 \
    gawk \
    build-essential \
    zlib1g-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    pandas \
    gffutils

# Install gffcompare (version 0.12.9)
RUN wget https://github.com/gpertea/gffcompare/releases/download/v0.12.9/gffcompare-0.12.9.Linux_x86_64.tar.gz && \
    tar -xzf gffcompare-0.12.9.Linux_x86_64.tar.gz && \
    mv gffcompare-0.12.9.Linux_x86_64/gffcompare /usr/local/bin/ && \
    mv gffcompare-0.12.9.Linux_x86_64/trmap /usr/local/bin/ && \
    rm -rf gffcompare-0.12.9.Linux_x86_64* && \
    chmod +x /usr/local/bin/gffcompare /usr/local/bin/trmap

# Install gffread (from gffread repository)
RUN git clone https://github.com/gpertea/gffread.git && \
    cd gffread && \
    make release && \
    mv gffread /usr/local/bin/ && \
    cd .. && \
    rm -rf gffread && \
    chmod +x /usr/local/bin/gffread

# Verify installations
RUN gffcompare --version && \
    gffread --version && \
    gawk --version && \
    python3 --version && \
    python3 -c "import pandas; print('pandas version:', pandas.__version__)" && \
    python3 -c "import gffutils; print('gffutils version:', gffutils.__version__)"

WORKDIR /usr/local/src

CMD ["bash"]

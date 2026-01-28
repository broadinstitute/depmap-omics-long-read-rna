FROM python:3.11.14

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/app \
    UV_VENV=/app/.venv

ENV PATH="$UV_VENV/bin:$PATH"

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    wget \
    ca-certificates \
    libgomp1 \
    gawk \
    build-essential \
    zlib1g-dev \
    git \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

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

WORKDIR /app

RUN pip install uv

COPY pyproject.toml uv.lock README.md ./
COPY combine_requantify_tools/__init__.py ./combine_requantify_tools/__init__.py
RUN uv sync --no-cache --no-group=dev

COPY combine_requantify_tools ./combine_requantify_tools

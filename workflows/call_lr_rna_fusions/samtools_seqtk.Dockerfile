FROM us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools:production

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    bc \
    git

RUN git clone https://github.com/lh3/seqtk.git
RUN cd ./seqtk && \
    make && \
    cd .. && \
    mv ./seqtk/seqtk /usr/local/bin/ && \
    rm -r ./seqtk

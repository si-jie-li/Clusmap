# clusmap with the motif-analysis toolchain (MEME-suite + bedtools + samtools).
# Lets users run the AME motif step locally without an HPC cluster.
#
#   docker build -t clusmap .
#   docker run --rm -v $PWD:/work -w /work clusmap clusmap-config show
#   docker run --rm -v $PWD:/work -w /work clusmap \
#       python -c "import clusmap as cm; cm.motif_pipeline('HM_ModGene.csv', run_mode='local')"
#
# Mount your genome FASTA / motif DBs as volumes and point clusmap-config at the
# in-container paths.

FROM python:3.11-slim

ARG MEME_VERSION=5.5.7

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential wget curl perl libxml2-dev libxslt1-dev zlib1g-dev \
        ghostscript bedtools samtools \
        libexpat1-dev libhtml-template-perl \
    && rm -rf /var/lib/apt/lists/*

# --- MEME-suite (provides `ame`) ------------------------------------------- #
RUN cd /tmp \
    && wget -q https://meme-suite.org/meme/meme-software/${MEME_VERSION}/meme-${MEME_VERSION}.tar.gz \
    && tar xzf meme-${MEME_VERSION}.tar.gz \
    && cd meme-${MEME_VERSION} \
    && ./configure --prefix=/opt/meme --enable-build-libxml2 --enable-build-libxslt \
    && make -j"$(nproc)" && make install \
    && rm -rf /tmp/meme-${MEME_VERSION}*
ENV PATH="/opt/meme/bin:/opt/meme/libexec/meme-${MEME_VERSION}:${PATH}"

# --- seqkit (fast FASTA toolkit, optional) --------------------------------- #
RUN cd /tmp \
    && wget -q https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz \
    && tar xzf seqkit_linux_amd64.tar.gz && mv seqkit /usr/local/bin/ \
    && rm -f seqkit_linux_amd64.tar.gz

# --- clusmap --------------------------------------------------------------- #
WORKDIR /opt/clusmap
COPY . /opt/clusmap
RUN pip install --no-cache-dir ".[all]"

WORKDIR /work
CMD ["python", "-c", "import clusmap; print('clusmap', clusmap.__version__, '+ MEME ame ready')"]

FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="DotMatch" \
      org.opencontainers.image.description="Deterministic known-target short-DNA assignment for CRISPR guides, barcodes, panels, and whitelist-style target sets." \
      org.opencontainers.image.source="https://github.com/dnncha/dotmatch" \
      org.opencontainers.image.url="https://github.com/dnncha/dotmatch" \
      org.opencontainers.image.documentation="https://github.com/dnncha/dotmatch#readme" \
      org.opencontainers.image.version="0.1.2" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.authors="Donncha O'Toole"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    python3 \
    python3-pip \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/dotmatch
COPY . .
RUN make clean && make && make shared

ENV DOTMATCH_LIB=/opt/dotmatch/libdotmatch.so
ENV PYTHONPATH=/opt/dotmatch/python
ENTRYPOINT ["/opt/dotmatch/dotmatch"]

FROM debian:bookworm-slim

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

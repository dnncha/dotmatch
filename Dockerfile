# Build with the packaging compiler flags (portable ISA), not host-tuned make.
FROM python:3.11-slim-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /build
COPY . .
RUN python -m pip wheel --no-deps --wheel-dir /wheels .

FROM python:3.11-slim-bookworm AS runtime
LABEL org.opencontainers.image.title="DotMatch" \
      org.opencontainers.image.description="Deterministic known-target short-DNA assignment for CRISPR guides, barcodes, panels, and whitelist-style target sets." \
      org.opencontainers.image.source="https://github.com/dnncha/dotmatch" \
      org.opencontainers.image.url="https://dotmatch.readthedocs.io/" \
      org.opencontainers.image.documentation="https://dotmatch.readthedocs.io/" \
      org.opencontainers.image.version="0.5.0" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.authors="Donncha O'Toole"
COPY --from=build /wheels /wheels
RUN python -m pip install --no-index --no-deps --no-cache-dir /wheels/*.whl \
    && rm -rf /wheels
ENV PYTHONDONTWRITEBYTECODE=1
WORKDIR /data
# The packaged Python dispatcher includes native counting AND assay/sensitivity/
# agent workflows. Running as the caller's UID is supported with --user.
ENTRYPOINT ["/usr/local/bin/dotmatch"]

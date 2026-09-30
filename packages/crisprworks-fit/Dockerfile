FROM python:3.12-slim AS builder
RUN apt-get update && apt-get install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /build
COPY . .
RUN python -m pip wheel --wheel-dir /wheels -c container-constraints.txt .

FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends libstdc++6 libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /wheels /wheels
RUN python -m pip install --no-index --find-links=/wheels crisprworks-fit \
    && rm -rf /wheels
WORKDIR /data
ENTRYPOINT ["crisprworks-fit"]
CMD ["--help"]

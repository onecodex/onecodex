# syntax=docker/dockerfile:1.7
FROM ghcr.io/astral-sh/uv:python3.13-bookworm-slim AS builder

ENV UV_LINK_MODE=copy \
    UV_NO_INSTALLER_METADATA=1

WORKDIR /app
COPY . .

RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --system ".[all]"

FROM python:3.13-slim-bookworm

RUN groupadd --system --gid 65532 nonroot \
 && useradd  --system --uid 65532 --gid 65532 --home /home/nonroot --create-home --shell /sbin/nologin nonroot

COPY --from=builder /usr/local/lib/python3.13/site-packages /usr/local/lib/python3.13/site-packages
COPY --from=builder /usr/local/bin/onecodex /usr/local/bin/onecodex

USER nonroot
WORKDIR /home/nonroot

ENTRYPOINT ["python3"]

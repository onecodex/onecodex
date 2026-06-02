# syntax=docker/dockerfile:1.7
FROM ghcr.io/astral-sh/uv:python3.11-bookworm-slim AS builder

ENV UV_LINK_MODE=copy \
    UV_NO_INSTALLER_METADATA=1

WORKDIR /app

COPY pyproject.toml README.md ./

RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --system ".[all]"

# Source layer: install just the project on top of the cached deps.
RUN rm -rf onecodex
COPY onecodex/ ./onecodex/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --system --no-deps --force-reinstall .

# Precompile .pyc up-front. This reduces start time by about ~260ms
RUN python -m compileall -q -j 0 /packages

FROM python:3.11-slim-bookworm

RUN groupadd --system --gid 65532 nonroot \
 && useradd  --system --uid 65532 --gid 65532 --home /home/nonroot --create-home --shell /sbin/nologin nonroot

COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages

USER nonroot
WORKDIR /home/nonroot

ENTRYPOINT ["python3"]

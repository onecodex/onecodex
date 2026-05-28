FROM python:3.11-slim-bookworm AS builder

WORKDIR /app
COPY pyproject.toml README.md ./
COPY onecodex/ ./onecodex/

RUN pip install --no-cache-dir --target=/packages ".[all]"


FROM gcr.io/distroless/python3-debian12:nonroot

COPY --from=builder /packages /packages
ENV PYTHONPATH=/packages

ENTRYPOINT ["python3"]

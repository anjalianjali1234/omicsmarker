FROM python:3.11.11-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/omicsmarker
COPY pyproject.toml .
COPY src src
RUN pip install --no-cache-dir '.[analysis]'

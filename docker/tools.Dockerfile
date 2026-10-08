# Small Python 3 tool image: converts the published xlsx into the canonical CSV read by both
# the legacy baseline and the port (python -m embd_data.sheet). openpyxl pinned.
FROM python:3.12.7-slim-bookworm
RUN pip install --no-cache-dir openpyxl==3.1.5
ENV PYTHONPATH=/repo/src PYTHONDONTWRITEBYTECODE=1
WORKDIR /repo

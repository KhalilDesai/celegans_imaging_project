# EMBD downstream analysis (Python 3 port): reproducible runtime image.
#
#   docker build -t embd:latest .                         # runtime image (code + locked dependencies)
#   docker build --target test -t embd:test .             # + tests, golden outputs and the Zenodo data
#
# Run (the input data is mounted read-only, outputs are written to a mounted directory):
#   docker run --rm -v "$PWD/data":/data:ro -v "$PWD/out":/out embd:latest
#
# Both base images are pinned by multi-architecture digest (works on linux/amd64 and linux/arm64).
# Python dependencies come from uv.lock (uv sync --frozen): the same versions on every build.

# --- base: Python 3.12 + locked dependencies -------------------------------------------------
FROM python:3.12.12-slim-bookworm@sha256:593bd06efe90efa80dc4eee3948be7c0fde4134606dd40d8dd8dbcade98e669c AS base
COPY --from=ghcr.io/astral-sh/uv:0.9.27@sha256:143b40f4ab56a780f43377604702107b5a35f83a4453daf1e4be691358718a6a /uv /usr/local/bin/uv

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PATH=/opt/venv/bin:$PATH \
    PYTHONPATH=/app/src \
    MPLBACKEND=Agg \
    # Determinism: single-threaded BLAS / numba (see docs/RANDOMNESS.md).
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    NUMBA_CACHE_DIR=/tmp/numba_cache MPLCONFIGDIR=/tmp/mpl

WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --no-install-project --no-dev

COPY src ./src

# --- runtime ----------------------------------------------------------------------------------
FROM base AS runtime
RUN useradd --create-home --uid 1000 embd && mkdir -p /out && chown embd /out
USER embd
VOLUME ["/data", "/out"]
ENTRYPOINT ["python", "-m", "embd.pipeline"]
CMD ["reproduce", "--sheet", "/data/All_genes_mySQL_table.xlsx", "--index", "/data/Index_EMBD_ID_Gene_names.xlsx", "--out", "/out"]

# --- test: adds pytest, the tests (incl. golden legacy outputs) and the Zenodo data -------------
FROM base AS test
RUN uv sync --frozen --no-install-project --extra test
COPY tests ./tests
COPY tools ./tools
COPY data ./data
COPY development ./development
COPY embryos_analyze ./embryos_analyze
RUN useradd --create-home --uid 1000 embd && chown -R embd /app
USER embd
ENTRYPOINT ["python", "-m", "pytest"]
CMD ["-q"]

# Default target is the runtime image.
FROM runtime

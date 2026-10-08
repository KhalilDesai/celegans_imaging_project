# Legacy baseline: runs the ORIGINAL Python 2.7 downstream code, unmodified, against the
# spreadsheet-backed MySQL emulation (src/embd_data). See legacy_baseline/README.md.
#
# Build (from repo root; amd64 to match the original x86-64 Windows workstation):
#   docker build --platform linux/amd64 -f docker/legacy.Dockerfile -t embd-legacy:py27 .
FROM python@sha256:d8fac68ebdc45b8d66d53f1ed6c1532da81109a8f5532a6ca0c951ed31107d70
# ^ python:2.7.18-buster (linux/amd64), pinned by digest

ENV PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLBACKEND=Agg \
    LC_ALL=C.UTF-8 LANG=C.UTF-8

# numpy first: matplotlib 1.4 builds from source against it.
RUN pip install numpy==1.13.3

# requirements.txt pins matplotlib==1.4.0rc2. Try it, then the nearest final 1.4.x release.
# The installed version is recorded in /legacy-versions.txt and in the run manifest.
ARG MPL_CANDIDATES="1.4.0rc2 1.4.0 1.4.3"
RUN for v in $MPL_CANDIDATES; do \
      pip install "matplotlib==$v" && echo "matplotlib==$v" > /mpl-version.txt && break; \
    done; test -s /mpl-version.txt

COPY legacy_baseline/requirements-legacy.txt /tmp/requirements-legacy.txt
RUN pip install -r /tmp/requirements-legacy.txt "$(cat /mpl-version.txt)" \
 && pip freeze > /legacy-versions.txt && python -V 2>> /legacy-versions.txt

RUN useradd --create-home --uid 1000 embd
USER embd
WORKDIR /work
# The repository is bind-mounted read-only at /repo and the output dir at /work at run time.
ENV EMBD_REPO=/repo EMBD_WORK=/work EMBD_SQLITE=/work/embd_mysql_emulation.sqlite
ENTRYPOINT ["python", "/repo/legacy_baseline/run_legacy.py"]

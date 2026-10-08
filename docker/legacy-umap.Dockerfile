# Legacy baseline, stage B: development/makeUMAPs_V2.py ran under Python 3.6 (its companion
# groupAssesment.make_UMAP_content_files says "needs python 3.6 to make the UMAP itself").
# The original umap-learn/numba versions are unknown (DECISIONS C1/C2, Q6); the pins below are a
# period-appropriate guess for April 2022 and are recorded with every run.
#   docker build --platform linux/amd64 -f docker/legacy-umap.Dockerfile -t embd-legacy:umap-py36 .
FROM python:3.6.15-slim-bullseye
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
RUN pip install numpy==1.19.5 scipy==1.5.4 scikit-learn==0.24.2 pandas==1.1.5 \
      matplotlib==3.3.4 seaborn==0.11.2 numba==0.53.1 llvmlite==0.36.0 pynndescent==0.5.6 umap-learn==0.5.3 \
 && pip freeze > /legacy-umap-versions.txt && python -V >> /legacy-umap-versions.txt
RUN useradd --create-home --uid 1000 embd
USER embd
WORKDIR /work
ENV EMBD_REPO=/repo EMBD_WORK=/work
ENTRYPOINT ["python", "/repo/legacy_baseline/run_umap_legacy.py"]

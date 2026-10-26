# Handoff: where the port stands (for a new Claude Code session or a new machine)

Read this first, then `docs/PORTING_PLAN.md`, `docs/DECISIONS_AND_ISSUES.md`,
`docs/NUMERICAL_DIFFERENCES.md` and `docs/PERFORMANCE.md`.

## The project in one paragraph
Python 2.7 → Python 3 port of the **downstream** analysis of Green et al. 2024 (*Cell*): a per-gene parameter
table → PAD similarity → best matches, reciprocal matches, auto-descriptions, UMAP. Image analysis is out of
scope (it raises `NotImplementedError`). The input is `data/All_genes_mySQL_table.xlsx` (Zenodo 7761843).
Scientific correctness comes first: the port must reproduce the original code's results.

## Status (2026-10-26)
| Phase | State |
|---|---|
| 0 Audit | done: `docs/PORTING_PLAN.md` |
| 1 Data layer + legacy baseline | done. The original Py2.7 code was run unmodified (Docker, `legacy_baseline/`); its outputs are frozen in `tests/golden/` (maintainer chose option a) |
| 2 Python 3 port | done: `src/embd/` (names mirror the originals), `src/embd_data/` (spreadsheet → SQLite MySQL emulation) |
| 3 Performance | done: full pipeline ≈ 145 s natively, outputs byte-identical to before the optimisations |
| 4 Reproducibility | done except the README (see below) |

## Phase 4: what exists and what's left
- Done: `uv.lock`; `Dockerfile` (base images pinned by digest, `runtime` and `test` targets, non-root user,
  data at `/data`, outputs at `/out`); `python -m embd.pipeline reproduce --sheet … --index … --out …`;
  `Makefile` (`make test`, `make reproduce`, `make docker-build`, `make docker-reproduce`, `make docker-test`, …);
  `apptainer.def`; `.github/workflows/ci.yml`; `docs/RANDOMNESS.md`; `tools/mem_profile.sh`.
- Memory (2026-10-26): every step peaks ≤ 0.46 GB natively (PERFORMANCE "Memory", #5). The earlier freeze was
  most likely the fixed 4 GB Docker VM on an 8 GB host. Still always cap containers (`make docker-reproduce`
  uses `--memory=3g`).
- Container verified on a 16 GB Apple-silicon Mac (Docker Desktop, 8 GB VM): full pipeline 236 s, peak 431 MiB;
  all 56 tests pass on its outputs. Native macOS and the image differ in the last bits (Accelerate vs OpenBLAS;
  NUMERICAL_DIFFERENCES N1 "Platform"); rankings identical. For byte-identical outputs, use the image.
- Not yet run: the Apptainer build (no Apptainer on the dev machines) and the CI workflow (runs on push; GitHub
  runners are amd64, so outputs may again differ in the last bits from arm64, within N1).
- **README: on hold.** The maintainer wants a guided walk-through of the code first.

## Working agreements with the maintainer
- Faithful first. Don't change scientific behaviour without approval. Log every finding/decision in
  `docs/DECISIONS_AND_ISSUES.md`. If a library upgrade changes numbers, stop and show the size.
- Before anything heavy (more than a minute, or lots of memory), check free memory and ask. Never run the full
  pipeline inside Docker on a machine with < 16 GB RAM without a `--memory` cap.
- Commit small and often, and push to `origin` (the maintainer's repo). Never touch `regreen/EMBD_public`.
- Open decisions: A18 (how to handle duplicated genes across the ~1,900-gene set), an optional 1–500-only mode,
  per-gene UMAPs for every query gene, and porting the non-published experiment code (tSNE, classifiers).

## Set up on a new machine
```bash
git clone https://github.com/KhalilDesai/celegans_imaging_project.git
cd celegans_imaging_project
curl -LsSf https://astral.sh/uv/install.sh | sh     # installs uv
uv sync --frozen --extra test                       # Python 3.12 + locked dependencies in .venv
.venv/bin/python -m pytest -q                       # fast tests (~40 s)
.venv/bin/python -m embd.pipeline reproduce --sheet data/All_genes_mySQL_table.xlsx \
    --index data/Index_EMBD_ID_Gene_names.xlsx --out out   # full pipeline (~2.5 min)
EMBD_PORT_RUN_DIR=out .venv/bin/python -m pytest -q -m slow   # compare every output with tests/golden
```
On macOS, if `import embd` fails, the venv's `.pth` file may have been flagged hidden:
`chflags nohidden .venv/lib/python3.12/site-packages/*.pth` (DECISIONS P2-8).

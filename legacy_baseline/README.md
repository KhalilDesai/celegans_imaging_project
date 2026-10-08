# Legacy baseline (original Python 2.7 code, unmodified)

This directory runs the **original** downstream code from `development/` and `embryos_analyze/`
under Python 2.7 to produce the reference ("golden") outputs that the Python 3 port is tested
against. **No original source file is edited.** The repository is mounted read-only into the container.

How the unmodified code is made runnable:

| Problem in the public repo | Fix (outside the original files) |
|---|---|
| Needs a MySQL server (`MySQLdb`, `genes`/`pad_rank`/`pair_lookup` tables) | `shims/MySQLdb.py` → `src/embd_data/mysql_emulation.py`: SQLite file built from the published sheet. Writes go to `build/legacy_work/embd_mysql_emulation.sqlite`. |
| `varLookup.MYSQL_DB_CREDENTIALS`, `db_utils_embryos.mySQLpath/login`, `fileLookup.SKIPPED` undefined | set by `run_legacy.py:bootstrap()` before use |
| Imports of image-analysis modules (`cv2`, `PIL`, `tifffile`, `pims`, `Embryos`, `emb_handler`, `avgCalculator`) and of modules missing from the repo (`dataValidation`) | stubs in `shims/`; calling them raises `NotImplementedError("image analysis not yet ported")` |
| `Z:/…` paths | the container's working dir contains a directory literally named `Z:` |
| Parameter set chosen by editing `get_params_use()` | `run_legacy.py` swaps it before `RNAiClass` is imported, one process per set |
| Hand-edited loop bounds in `update_top_hits_PAD`, `update_column_reciprocal_matches` | verbatim copies of the bodies in `run_legacy.py` with loop = EMBD 1–503 |

Every such deviation is listed in `docs/DECISIONS_AND_ISSUES.md`, section L.

## Run

```bash
docker build --platform linux/amd64 -f docker/legacy.Dockerfile -t embd-legacy:py27 .
docker build -f docker/tools.Dockerfile -t embd-tools:py312 .
docker run --rm -v "$PWD":/repo:ro -v "$PWD/build/canonical":/out embd-tools:py312 \
    python -m embd_data.sheet data/All_genes_mySQL_table.xlsx data/Index_EMBD_ID_Gene_names.xlsx /out
legacy_baseline/run.sh init-db
for s in opt gls ms all; do legacy_baseline/run.sh phenodist $s; done   # slow (CSI is ~O(n^4))
legacy_baseline/run.sh bestmatches
legacy_baseline/run.sh padrank
legacy_baseline/run.sh derived
for s in opt gls ms all; do legacy_baseline/run.sh matrices $s; legacy_baseline/run.sh vectors $s; done
```

Outputs land in `build/legacy_work/Z:/`. Step timings go to `Z:/timings_legacy.jsonl`.
`compare_reference.py` (Python 3) compares the rankings with the published `all_best_matches_05192022.xlsx`.

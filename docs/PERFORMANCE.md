# Performance

Machine: Apple-silicon laptop (8 GB RAM, of which ~4 GB was held by the Docker VM), Python 3.12, numpy 2.5.
Wall-clock seconds as logged by `embd.pipeline` (`timings_port.jsonl`). phenoDist uses the original
`multiprocessing.Pool(6)`. Every optimisation is its own commit and must keep the full regression suite
(`pytest` and `pytest -m slow`) passing.

## Starting point

| Step | Legacy (Py2.7, amd64 emulation) | Port, Phase 2 (before Phase 3) |
|---|---|---|
| phenoDist (all pairs: dist, PAD, CSI), per parameter set | 2,630–2,890 s | 66–77 s |
| best matches: `output_all_best_matches_to_csv_combined` | 24,032 s (≈6.7 h; 40 h on a rerun, see X8) | 107 s |
| best matches: `output_all_best_matches_to_one_csv_combined` | 27,328 s | 94 s |
| PAD_RANK | 2,800 s | 0.8 s |
| derived columns (descriptions, WT-like, top hits, reciprocal) | ≈270 s | 13 s |
| matrices for UMAP (PAD + distance, 503²), per parameter set | ≈4,100–4,600 s | 46–62 s |
| **whole pipeline** | **≈ 2–3 days** | **≈ 14 min** |

The Phase-2 gains came from output-neutral caching (DECISIONS P2-4) and native execution.

## Profile (Phase 3, before any change)

`tools/profile_steps.py`, single process, cProfile, optimal parameter set, extrapolated to 503 genes:

| Step | Single-core cost | Where the time goes |
|---|---|---|
| all-pairs PAD (`getPADS`) | ≈ 78 s | `getNDDistance` is called **5× per pair** (1.27 M calls). 3 of the 5 recompute a gene's distance to control (`getDist2Zero`, 759 k calls, 56 s), which never changes. Builtin `sum()` over a numpy bool array costs 16 s. `getOverlapDims` rebuilds four filtered lists per pair (8 s). |
| CSI (`getCSIS` → `getNBelow`) | ≈ 94 s | For every (i, j): a linear scan over all 503 genes to find gene j by name, then `r2 in rcounted` membership tests against a Python **list**. |
| matrices (`make_UMAP_content_files`) | ≈ 53 s | Recomputes every PAD and distance through the same per-pair methods (same `getDist2Zero` waste). |
| best matches | ≈ 15 s per 20 seeds (incl. a one-off 5 s CSV load) | The label index built for `PADFilter` (P2-4) is **rebuilt on every call**, because the original re-slices `readCSV(...)[3:]` each iteration, so the list object (the cache key) is new each time (8.4 s). `getGeneCommonName` re-reads `gene_key.csv` on every call (≈10 per seed). |

## Optimisations

| # | Change | Output | Before | After |
|---|---|---|---|---|
| 1 | Cache each gene's distance to control (`getDist2Zero`), cleared by `setNDPosition`; `not vOverlap.any()` instead of builtin `sum()` over a numpy array | **byte-identical** (phenoDist CSV, PAD/Dist matrices) | phenoDist opt 70.8 s; matrices opt 45.8 s | phenoDist opt **38.0 s**; matrices opt **8.2 s** |
| 2 | CSI: same counts computed with array operations per seed (distance matrix reordered by the seed's ranking, strict lower-triangle row sums), instead of a name scan over all genes plus list membership per pair; original kept as `getCSI_original`/`getCSIS` | **byte-identical** (phenoDist CSV); equality also tested on synthetic data full of ties/inf | phenoDist opt 38.0 s | phenoDist opt **8.8 s** |
| 3 | Best matches: `PADFilter` row index keyed on the shared row objects (the original's per-iteration `[3:]` re-slice defeated the list-keyed cache); gene-name lookups read `gene_key.csv` once per file version | **byte-identical** (`all_best_matches`, `hits_matrix_v2`) | combined 107 s; one-csv 94 s | combined **14.0 s**; one-csv **10.2 s** |

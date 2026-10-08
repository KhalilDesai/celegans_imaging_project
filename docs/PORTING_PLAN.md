# EMBD Downstream Analysis — Porting Plan (Phase 0 audit)

Source: fork of [regreen/EMBD_public](https://github.com/regreen/EMBD_public) (Python 2.7), accompanying
Green et al. 2024, *Cell* 187:3141, "Automated profiling of gene function during embryonic development"
(DOI [10.1016/j.cell.2024.04.012](https://doi.org/10.1016/j.cell.2024.04.012), PMC11166207).

Companion log: [DECISIONS_AND_ISSUES.md](DECISIONS_AND_ISSUES.md). It holds every factual finding,
scientific concern, and design decision, each with an ID such as `B8`. This plan cites those IDs instead of repeating them.

**Porting policy** (maintainer, 2026-10-08): port faithfully and bug-for-bug first, in its own commit(s). Anything
that looks wrong is logged, not fixed. Performance work comes only after the faithful port passes its
regression tests, and it must keep results identical within documented tolerance.

---

## 1. Scope

### 1.1 In scope (downstream: per-gene parameter table → PAD → rankings → groups/UMAP → figures)

| Module | Role | Live entry points |
|---|---|---|
| `embryos_analyze/varLookup.py` | `PARAM_NAMES` (138), constants, `printLog` | — |
| `development/fileLookup.py` | Input/output paths (all under `Z:/`) | — |
| `development/params_use_assign.py` | Parameter subsets: `_group_rank_opt` (optimal 17+24), `_all` (58+50), `_gls_params`, `_ms_params` | `get_params_use()` |
| `development/embdFunc.py` | `getNDDistance`, `getNDAvgVector`, `split_params`, gene-name lookups | — |
| `development/RNAiClass.py` (subset) | `RNAiClass.__init__`, `load_params_data`, `setOrigin` (hardcoded constants), `setNDPosition`, `getDistance`, `getDist2Zero`, `getPAD`, `get_cosine_dist`, `get_norm_dist`, `getSimParams`, `getDiffParams`, `getOverlapDims`, `getDims`, `ControlClass.__init__` | — |
| `development/phenoDistCalc.py` | All-pairs distance/PAD/CSI → `phenoDist_rank_MMDDYYYY.csv` | `run_phenodistcalc_populate()` |
| `development/PADassist.py` (subset) | `PADFilter`, `get_data*`, `update_pad_rank` (→ PAD_RANK), best-match CSV writers, `makeNetwork`, `getCommonPairs`, `getPairsParams` | `output_all_best_matches_to_csv_combined`, `output_all_best_matches_to_one_csv_combined`, `update_pad_rank` |
| `development/db_utils_genes.py` (read paths + derived writers) | `get_id_from_embd`, `get_reciprocal_matches`, `update_top_hits_PAD`, `update_column_reciprocal_matches`, `update_WT_like_or_phenotype`, `update_automated_description` | the `update_*` functions (they write to the DB, which will be redirected to local files) |
| `development/groupAssesment.py` (subset) | Manual groups, centroids, closest group / `auto_description`, PAD/dist matrices, tSNE, UMAP input files, rank metrics | `get_all_closest_groups`, `make_UMAP_content_files`, `make_tSNE_all_data` |
| `development/makeUMAPs_V2.py` | UMAP figures (originally run under py3.6, see C1) | script |
| `development/cumulativeDistribs.py` (subset) | `getAllPCC`, `load_pcc`, `get_pca`, `get_lda`, `plotCum`, `plot_scikit_lda` | — |
| `development/make_figures.py` (subset) | `plot_pretty_graphs`, `show_500gene_heatmap`, `show_manual_heatmap_by_gene`, `plot_tIntVsPtArrest`, `plot_maxIntGLSVsPtArrest`, `sigmoidGLS/MS` | — |
| `embryos_analyze/myFigure.py` | matplotlib wrapper (image-free) | — |
| `embryos_analyze/myFunc.py` (subset) | `readCSV`, `sort_nicely`, `mkdir_p`, `find_consec_inds` | — |
| `embryos_analyze/myMath.py` (subset) | `getSig`, `getSubIndex`, `nan_helper` | — |
| `embryos_analyze/db_utils_embryos.py` (subset) | `initialize_db`, `get_row`, `get_columns`, `update_row`, `insert_row`, which are used generically for the `genes`/`pad_rank`/`pair_lookup` tables | replaced by the data layer (§4) |

### 1.2 Out of scope (left byte-for-byte unchanged)

- `embryocrop/` (entire, including the PyQt GUI).
- `embryos_analyze/`: `Embryos.py`, `MSImageClass.py`, `timeAlign.py`, `moveDetect.py`, `avgCalculator.py`,
  `emb_handler.py`, `db_utils_embryos.py` write paths for the `embryos` table, image functions of
  `myFunc`/`myMath`, `UnitTests_Embryos/`.
- `RNAiClass`: `refresh_params_data`, `refresh_embs`, `getEmbryos`, `splitByMovement`, `calcAverageParams`,
  `asignUse`, `load_from_pickle`, `getPrev`, `getEmbsUsed`, `getRepEmb`, all `plot*`/`showMaxProj*`/`makeTalignMovies`,
  `fix_head_params_RNAi`, `ControlClass.defineFolders/getEmbryos/asignUse`, `get_control_origin`, `write_controlAvg_mysql`.
- `PADassist`: movie functions (`show_best_match_movies`, `make_movies`, `convert_tifstack2mp4`), Yanai lookups (CSV
  not provided), `findCliques`/`getInteractions` (need missing `dataValidation`, interactome CSV).
- `make_figures`: every function that loads Embryo objects or images (L17–881, L1097–2317, see U2).
- `cumulativeDistribs.show_param_dist/get_pca_control*`, `paramDistribs` (all per-embryo; see P7/U2).
- `development/misc_programs/`, `UT_bench_mark_emb.py`, `UnitTests/UT_rnai_save_load.py`.

**Layout proposal** (for approval): the ported package lives in `src/embd/`. All original files stay where they
are, unmodified, and become the "legacy" reference. Nothing is moved, so no deletions are needed. A `legacy/README.md`
explains that `embryocrop/`, `embryos_analyze/` (image parts), and the image functions in `development/` are not ported yet.
Alternative: `git mv` the out-of-scope code under `legacy/`. That is still byte-identical, but it changes paths
that the legacy container would then need to be told about. I recommend leaving files in place.

### 1.3 Import graph (module level; ★ = missing from repo, ⊘ = out of scope / image)

```
varLookup ── (none)
fileLookup ── varLookup
params_use_assign ── numpy
myFigure ── matplotlib, scipy.stats
myFunc ── ⊘cv2, ⊘tifffile, ⊘PIL, scipy.ndimage, myFigure
myMath ── scipy, matplotlib, lmfit, myFigure, myFunc
embdFunc ── myFunc, myMath, varLookup, fileLookup, numpy, pandas, git(GitPython)
db_utils_embryos ── MySQLdb, varLookup(MYSQL_DB_CREDENTIALS — not defined: commented out in varLookup!)
db_utils_genes ── pandas, db_utils_embryos, fileLookup, embdFunc
RNAiClass ── varLookup, fileLookup, ⊘emb_handler(→⊘Embryos→★SkeletonAnalysis, ⊘MSImageClass, ⊘timeAlign, ⊘moveDetect, skimage, cv2),
             myFigure, myMath, myFunc, ⊘cv2, ⊘PIL, db_utils_embryos, db_utils_genes, embdFunc, params_use_assign
PADassist ── myFunc, embdFunc, fileLookup, sklearn.metrics, myFigure, varLookup, RNAiClass, networkx, myMath,
             ★dataValidation, ⊘Embryos, ⊘cv2, ⊘PIL, db_utils_embryos, ⊘avgCalculator
phenoDistCalc ── varLookup, params_use_assign, RNAiClass, PADassist, lmfit, embdFunc
groupAssesment ── myFunc, varLookup, RNAiClass, PADassist, myFigure, embdFunc, ★dataValidation, sklearn, params_use_assign
cumulativeDistribs ── RNAiClass, myFigure, scipy, sklearn, ⊘avgCalculator, groupAssesment, myFunc, params_use_assign,
                      db_utils_embryos, ⊘emb_handler
makeUMAPs_V2 ── numpy, matplotlib, seaborn, umap, pandas   (standalone; ran under py3.6)
```

**As shipped, none of `RNAiClass`, `PADassist`, `phenoDistCalc`, `groupAssesment` can be imported from this repo.**
Reasons: the missing `SkeletonAnalysis`/`dataValidation` modules (C6), and `varLookup` never defines
`MYSQL_DB_CREDENTIALS` (it is commented out, so `db_utils_embryos` fails with an ImportError). The legacy baseline (§5)
fixes this **without editing any original file** by putting shim/stub modules first on `PYTHONPATH`.

### 1.4 Every point where downstream code reaches into image analysis or MySQL

| Call site | Reaches into | Port treatment |
|---|---|---|
| `RNAiClass.__init__` → `db_utils_embryos.initialize_db`, `db_utils_genes.get_id_from_embd` | MySQL `genes` | data layer (sheet) |
| `RNAiClass.load_params_data` → `get_row(genes, PARAM cols + prev_GLS, prev_MS, version)` | MySQL `genes` | data layer; supplies `version` (A4) |
| `load_params_data` → `load_from_pickle` → `refresh_params_data` (when `version is None` or `check_version`) | pickles, embryos, images, git | `NotImplementedError("image analysis not yet ported")`; never triggered when `version` is supplied |
| `RNAiClass.save_params_data` | MySQL `UPDATE genes` | local-file writer (never called on the faithful path) |
| `phenoDistCalc.write_phenodistinfo_to_mySQL` | MySQL `pair_lookup` | local file (`pair_lookup.csv`) |
| `PADassist.get_data_from_mysql` | MySQL `pair_lookup` | data layer reads `pair_lookup.csv` |
| `PADassist.update_pad_rank` | MySQL `pad_rank` (FLOAT pads) | local `pad_rank.csv` + float32 emulation (B8) |
| `db_utils_genes.get_reciprocal_matches`, `update_top_hits_PAD` | MySQL `PAD_RANK` | read local `pad_rank` |
| `db_utils_genes.update_column_reciprocal_matches`, `update_WT_like_or_phenotype` | MySQL `genes` (movement, WT_like, auto_description) | read sheet; write local `genes_derived.csv` |
| `groupAssesment.find_closest_groups` | MySQL `genes` read/write `auto_description` | read sheet; write local |
| `groupAssesment.retrieve_and_map_man_group_autodata` → `get_sql_data_for_man_group_conditions` | MySQL `genes` (`man_group`, `man_*`) | sheet (columns present) |
| `cumulativeDistribs.show_param_dist`, `get_control_origin`, `get_all_params` | MySQL `embryos` table (per-embryo) | out of scope → `NotImplementedError` |
| `embdFunc.get_rnai_revision_version` (`check_version=True`) | git repo of the analysis code | not called (default False) |
| `PADassist`, `RNAiClass`, `make_figures` movie/figure code | images via cv2/PIL/Embryos | `NotImplementedError` |

---

## 2. Pipeline: entry points and run order

The original pipeline was run by hand, by editing `__main__`/`initialize()` blocks and commented-out alternatives.
Reconstructed order for the published May 2022 outputs:

```
(0) genes table (sheet)                         [boundary — see §4 interface]
     │  RNAiClass(i).load_params_data()          for i in 1..503        (B10: 10 extra IDs not in loops)
     │  setOrigin(): hardcoded controlOrigin / rnaiStd                   (R1)
     │  setNDPosition(): posND over active parameter set                 (B2, B3)
(1) phenoDistCalc.run_phenodistcalc_populate()
     │  getPADS (pool of 6) → per i: rows (r1, r2, dist, dist2WT, dims, PAD) sorted by dist desc (B7)
     │  getCSI → CSI per pair  (O(n³–n⁴); see §7)
     └→ phenoDist_rank_04142022.csv                      [params = _group_rank_opt]
        phenoDist_rank_GLS_ONLY05102022.csv              [params = _gls_params]   (edit get_params_use)
        phenoDist_rank_MS_ONLY05102022.csv               [params = _ms_params]
        phenoDist_rank_ALL_PARAMS05102022.csv            [params = _all]
(2) PADassist.output_all_best_matches_to_csv_combined()  → all_best_matches_MMDDYYYY.csv  (top 10 per set)
    PADassist.output_all_best_matches_to_one_csv_combined() → hits_matrix_v2_MMDDYYYY.csv (PAD>0.5, 5-char PADs)
    PADassist.update_pad_rank(data_opt)                   → pad_rank table (top 50, float32)      (B8)
(3) db_utils_genes.update_WT_like_or_phenotype           → genes.WT_like_or_phenotype (movement + auto_description)
    db_utils_genes.update_top_hits_PAD                   → genes.top_hits_PAD (top 10, PAD>0.5)
    db_utils_genes.update_column_reciprocal_matches      → genes.reciprocal_matches (B9, P4)
    db_utils_genes.update_automated_description → groupAssesment.get_all_closest_groups (EMBD 37–503)  (G2–G6)
(4) groupAssesment.make_UMAP_content_files()             → PAD/Dist matrices .npy (503×503) + group lists
    makeUMAPs_V2.py (py3.6, umap-learn)                  → UMAP_EMBD_042222_PAD_*.svg               (C1, C2, P6)
    groupAssesment.make_tSNE_all_data()                  → X_embedded_500genes_71076_*.csv + figure (G8)
(5) figures: closest-group bar charts, clustermap/heatmaps, PCC/PCA/LDA (subset in §1.1)
```

Note the circular dependency in step 3: `WT_like_or_phenotype` reads `auto_description`, which step 3's last item
writes. The published state depends on the order in which those functions were run (and re-run). OPEN: confirm with the lab. The faithful pipeline will
run them in the order above. The reference CSV `EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv`
will show which order matches.

**Reference-comparison mapping:**
- `all_best_matches_05192022.xlsx` ↔ step 2's `output_all_best_matches_to_csv_combined` (four parameter sets; the brief
  mentions three, but the code writes optimal, GLS, MS, and all-params).
- `EMBD0001_0504_…_052422.csv` ↔ step 3 (top hits, reciprocal matches, auto descriptions).

### 2.1 PAD_RANK (asked in the brief)

PAD_RANK is **a derived cache, not an input**. `PADassist.update_pad_rank` fills it from the optimal-parameter
`phenoDist_rank_04142022.csv` via `PADFilter(i, data, 51)` for i in 1..503. It stores `rank1..50` (int) and `pad1..50`
(MySQL `FLOAT`, so float32). It is read only by `get_reciprocal_matches` (seed side: `pad > 0.55` on the
float32 value; reverse side: rank membership only) and `update_top_hits_PAD` (`pad > 0.5`).
**The published reciprocal matches and top-hit columns depend on it**, but it is fully reproducible from the PAD matrix.
The port regenerates it in memory. To be faithful, it rounds PADs to float32 before threshold tests and keeps the
default-0 fill. Phase 1 will report how many seeds have a PAD within float32 rounding of 0.55 or 0.5, i.e. whether
the float32 emulation actually matters.

---

## 3. Dependencies actually used by in-scope code

| requirements.txt entry | In-scope use | Downstream container |
|---|---|---|
| numpy 1.13.3 | everywhere | ✔ |
| scipy 1.0.0 | stats (KS, pearsonr), curve_fit, ndimage (myFunc import only) | ✔ |
| pandas 0.21.1 | sheet/CSV I/O, groupAssesment | ✔ |
| scikit-learn 0.18.1 | TSNE, PCA, LDA, train_test_split, LogisticRegression, metrics | ✔ |
| matplotlib 1.4.0rc2 | figures | ✔ (1.4.0rc2 likely uninstallable from PyPI; legacy uses closest 1.4.x — VERIFY) |
| seaborn 0.8.1 | clustermap/heatmap, UMAP style | ✔ |
| networkx 1.11 | `PADassist` module import, `makeNetwork` | ✔ |
| lmfit 0.9.1 | `phenoDistCalc` module import (`optimizePAD` only), `myMath` import | ✔ (or lazy-import) |
| GitPython 2.1.8 | `embdFunc` module import (`get_rnai_revision_version` only) | lazy-import; not required |
| openpyxl 2.4.2 | none originally; **new**: reading the sheet | ✔ |
| **umap-learn** | `makeUMAPs_V2` | ✔ **not in requirements.txt** (C1) |
| mysqlclient / MySQL-python / mysql-connector | DB | ✗ replaced by data layer |
| pyensembl | `db_utils_genes.get_ensemble_ref` (lazy, not on pipeline) | ✗ |
| Pillow, tifffile, PIMS | image/movie functions only | ✗ |
| opencv-python, scikit-image, mahotas | image only | ✗ **confirmed not needed** |
| PyQt4, Theano, Lasagne, scikit-neuralnetwork, py2exe, pywin32/pypiwin32, PyInstaller, pefile, macholib, altgraph | none (embryocrop/packaging or unused) | ✗ **confirmed not needed** |
| basemap, Django, dask, intermine, vboxapi, PGen, dis3, image, sklearn-contrib-lightning, lmfit deps etc. | none | ✗ |

The legacy container (§5) needs the same set (pinned versions) plus shim modules. It does not need cv2, skimage, or MySQL,
because the out-of-scope modules are stubbed on `PYTHONPATH`. That makes the legacy baseline cheaper and avoids
building opencv 3.1 for Python 2.7.

---

## 4. The boundary interface: per-gene parameter table

This becomes `docs/PARAMETER_TABLE.md` in Phase 1. It is summarized here for review.

- **One row per RNAi condition**, keyed by `embd_id` (int; 0 = control). In the 2022 data, the conditions are 1–503
  plus 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905.
- **Parameter columns:** `{p}_GLS` and `{p}_MS` for `p` in `varLookup.PARAM_NAMES` (138 names; note that names such as
  `CoM0Rtd_GLS` already carry a strain suffix, giving columns like `CoM0Rtd_GLS_GLS`). Columns may be absent and are then
  treated as NaN (A5). The minimum required set is the union of the four parameter subsets (58 GLS + 50 MS names, listed in
  `params_use_assign`).
- **Semantics:** each value is the **median across the condition's usable embryos** of that strain. "Usable" means
  non-moving embryos if more than one embryo failed to move, otherwise all non-dead embryos (P2). Only embryos whose
  per-parameter fit passed `Embryo.checkGoodParam` contribute. Units are those of the feature-curve parameterization (time in
  20-min frames scaled by `tScale`, intensities in camera counts, distances/moments normalized by embryo radius/half-length,
  counts in nuclei, fractions dimensionless). Per-parameter units go in PARAMETER_TABLE.md.
- **NaN convention:** NaN/empty/`'NULL'` = parameter unavailable for that condition (e.g. marker not expressed). It
  is **not** zero. 0 is a legitimate value (A7). The distance metric handles NaN via overlap and penalty (B1).
- **Precision:** float64 in the table. The 2022 values were stored as MySQL float32 (A6). A producer should write
  float64 and document whether the values were rounded.
- **Metadata columns used downstream:** `movement_GLS`, `movement_MS` (fraction moving, but compared to 0/1 in code),
  `man_group` and `man_*` (classifier experiments), `gene_name`/`ce_id` (labels), `auto_description`/`WT_like_or_phenotype`
  (derived and written by the pipeline, but read back in step 3).
- **Data-layer extras** (not science): `version` (constant string, e.g. `"sheet-zenodo-7761843"`) and `prev_GLS`/`prev_MS`
  (NaN). These exist so that `load_params_data` never falls back to pickles or embryos (A4).

---

## 5. Phase 1 plan (data layer + legacy baseline)

1. **Data layer** (`src/embd_data/`, py2/py3-compatible pure Python + pandas/openpyxl). It provides drop-in replacements
   with **identical signatures** for `db_utils_embryos.initialize_db/get_row/get_columns/update_row/insert_row`, plus
   `db_utils_genes.get_id_from_embd` and the PAD_RANK/pair_lookup reads. They are backed by the sheet (read-only) and by local
   CSV files for writes. `initialize_db()` returns a fake `(conn, cursor)` that understands only the handful of SQL
   statement shapes actually used. Any other SQL raises loudly.
   - Cast all parameter columns to float64, and map `'NULL'`/empty to NaN.
   - Option to round parameter values to float32 to emulate the MySQL `FLOAT` columns (A6). The default is decided after
     the Phase 1 comparison.
   - `gene_key.csv` replacement built from the sheet or the Zenodo index (A9).
2. **Legacy container** `docker/legacy.Dockerfile`:
   - Stage A: `python:2.7` (pinned by digest), with pinned requirements for the in-scope set (§3).
   - On `PYTHONPATH` before the original dirs: the data-layer shims plus stubs for `cv2`, `PIL`, `tifffile`, `MySQLdb`,
     `skimage`, `SkeletonAnalysis`, `dataValidation`, `emb_handler`, `Embryos`, `avgCalculator`, and a
     `varLookup`-compatible credentials name. Each stub raises `NotImplementedError` on use.
   - `FOLDER_IN='Z:/'` is redirected by mounting the data/output directory at a path the code can open. If needed, a
     shim `fileLookup`/`varLookup` overlay changes **only** path constants. That is logged.
   - **Original `.py` files are not edited.** Drivers live in `legacy_runner/` and call the original functions in the
     order of §2, swapping `params_use_assign.get_params_use` via monkeypatch for the GLS/MS/all runs.
   - Stage B: `python:3.6` + umap-learn (version to be chosen near April 2022, e.g. 0.5.2/0.5.3, and **logged as a
     guess**) for `makeUMAPs_V2.py`.
3. **First experiment:** legacy PAD for all genes, compared with `all_best_matches_05192022.xlsx`. Report per
   parameter set (optimal, GLS, MS, all): exact rank agreement (top-10 set and order), Spearman/Kendall on the
   overlapping lists, and max/median |ΔPAD|. Run it with float64 sheet values and again with float32-rounded values. **Stop and
   report before freezing goldens.**
4. **Goldens** (after approval): `tests/golden/` with PAD and distance matrices (`.npy`, float64), the four phenoDist CSVs,
   rankings, PAD_RANK, top hits, reciprocal matches, closest groups and descriptions, UMAP/tSNE coordinates, and
   `MANIFEST.sha256`.

## 6. Phase 2 plan (py3 port)

- Target: Python 3.12, which every in-scope dependency supports (numpy 2.x, scipy, pandas 2.x, sklearn 1.5+,
  matplotlib 3.9+, seaborn 0.13, networkx 3, umap-learn 0.5.x + numba, lmfit 1.3). Python 3.13 is possible if numba and
  umap-learn wheels allow it at port time. VERIFY.
- Order: leaves first, i.e. `varLookup` → `params_use_assign` → utilities (`readCSV`, `getSig`, …) → `embdFunc` →
  `RNAiClass` (subset) → `phenoDistCalc` → `PADassist` → `db_utils_genes` (derived writers) → `groupAssesment` →
  `makeUMAPs_V2` → figures.
- Each module keeps its original public names. A thin package API is added later, not in the faithful commit.

### 6.1 Python 2→3 hazards in in-scope code

| Hazard | Where | Effect on results | Treatment |
|---|---|---|---|
| `print x` statements, `except X, e` | varLookup:118, embdFunc:86, PADassist (many), RNAiClass:1203+, phenoDistCalc:165,195, groupAssesment:111,1134,1149,1186,1430, cumulativeDistribs | syntax only | mechanical |
| `basestring` | db_utils_genes:535,636 | none | `str` |
| Integer division | **none on the PAD path** (all use `1.*` or float operands). make_figures `rows = len/cols` (out of scope); myFunc test-image fns (out of scope) | none in scope | audited; `//` only where py2 semantics were int (none in scope) |
| Dict ordering → vector order | `posND` uses `PARAM_NAMES` list order; `paramsUseGLS/MS` and `controlOrigin` dict order only used for membership | **none** (B11) | keep lists; tests assert parameter-vector order |
| Dict ordering → output order | CSV/DataFrame rows built from lists; `insert_row/update_row` SQL column order (irrelevant) | none | — |
| `map/zip/filter` as lists | groupAssesment:1959 (commented); make_figures:2307 (out of scope) | none in scope | `list(...)` |
| Float → str formatting | `phenoDistCalc.write_phenodistinfo_to_CSV` `'{}'.format(np.float64)`; PADassist best-match writers; `str(val)[0:5]` in hits_matrix_v2 | **yes**: CSV text and 12-digit rounding feed the ranking (B6). `str(val)[0:5]` truncation differs when py3 repr differs | explicit `'%.12g'` to reproduce py2 `str()`; test byte-equality of CSVs |
| Sorting with NaN / mixed types | `PADFilter` argsort on floats (PAD never NaN, B13); `load_pcc` sort with NaN; `getPADS` list sort by dist (inf allowed, no NaN possible?) | only in ties | `kind='stable'` is **not** a faithful choice; keep default and test for ties (B5) |
| Ties in rankings | `PADFilter` (argsort quicksort), `sort_values` in closest groups, `getPADS` stable `list.sort` | possibly | Phase 1 counts ties; if none in top 51 then irrelevant |
| `string` module | only `ascii_letters` imports (unused) | none | drop the import |
| str/bytes | `readCSV` text mode; csv `'wb'` in cumulativeDistribs; py2 pickles (out of scope); `gene_key.csv` latin-1 | encoding only | `open(..., newline='', encoding='latin-1')` where read |
| `np.random.seed(071076)` octal | groupAssesment:370 (commented) | none | — |
| Bare `raise` | embdFunc:29 | crash either way | `ValueError` (R5) |
| `is not 'dev delay'` | groupAssesment:1471 | none (prior branch) | `!=` |

### 6.2 Library-version hazards

| Library | Change | Affects | Treatment |
|---|---|---|---|
| numpy ≥1.24 | ragged `np.array([...])` raises | `RNAiClass.pNamesUse` (B12), `myFigure` L212 | `dtype=object` |
| numpy ≥1.24 | `np.float`, `np.int` removed | phenoDistCalc:49, groupAssesment:964 | builtins |
| numpy 1.25/2.x | SIMD sort kernels (x86-simd-sort AVX-512/AVX2; Highway on ARM) change **tie order** of default `argsort` | PADFilter rankings | test for ties; document |
| numpy 2.x | scalar repr `np.float64(0.5)`; NEP 50 promotion | CSV text via `str()`; float32×python-float promotion in float32 emulation | explicit formatting; explicit dtypes |
| numpy | `nanmean` of all-NaN → warning + NaN (relied on) | centroids | keep, silence warning |
| pandas 2.x | `read_csv(f, ',')` positional sep removed | RNAiClass:1787, db_utils_genes:101,846 | `sep=','` |
| pandas 2.x | `idxmax` on all-NaN warns or raises; downcasting FutureWarnings; copy-on-write | closest-group description | replicate the old NaN result explicitly |
| pandas | default `sort_values` kind quicksort (unstable) | closest group ties | keep and test |
| scipy ≥1.3 | `ks_2samp` default exact for small n | KS experiments (not in published outputs) | pass `method='asymp'` to reproduce |
| scipy | `scipy.optimize.minpack`, `scipy.ndimage.filters`, `scipy.stats.stats` private namespaces | imports | public namespaces |
| sklearn 0.19/1.2/1.5 | TSNE: early_exaggeration 4→12, init random→pca, learning_rate 200→auto, n_iter→max_iter | tSNE coordinates | pass the 0.18.1 defaults explicitly; coordinates compared by structure (trustworthiness/Procrustes), not bitwise |
| sklearn 0.22 | LogisticRegression solver liblinear→lbfgs | `log_regression` demo (not published) | pass `solver='liblinear'` |
| sklearn ≥0.24 | LDA `n_components` > classes−1 raises | `get_lda_by_groups` (already broken, G7) | port as-is, document |
| sklearn | `train_test_split` with global RNG seeded 1 | `getRNAiGroupsPheno` split (not on published path) | test that the split equals the legacy split |
| umap-learn 0.4/0.5 | defaults and parallel nondeterminism; random_state=None uses global RNG | UMAP coords | C2: faithful = same call; reproducible mode needs approval |
| matplotlib 2.x/3.x | `normed`, `set_axis_bgcolor`, `'on'/'off'` ticks, `tick.label`, `edgecolor=''`, `transparant=` typo | figure rendering only | per-call-site fixes, logged (U4) |
| seaborn 0.13 | `sns.set` → alias of `set_theme`; heatmap/clustermap compatible | figure styling | keep |
| networkx 2/3 | `cliques_containing_node` removed in 3.x | `findCliques` (out of scope: needs dataValidation) | — |
| lmfit 1.x | Parameters API stable; `minimize(..., options=)` | `optimizePAD` (not on published path) | port as-is |

## 7. Hardcoded paths in in-scope code

All are rooted at `FOLDER_IN = 'Z:/'` (varLookup:13) or literal `'Z:/…'`/`'Z:\\…'`. Count by file (non-comment lines):
PADassist 28, groupAssesment 18, make_figures 18, makeUMAPs_V2 9, phenoDistCalc 3, cumulativeDistribs 2, embdFunc 1,
RNAiClass 1, varLookup 1. The faithful port routes every one of them through a single configurable root
(`EMBD_DATA_DIR` / `EMBD_OUT_DIR`) while keeping the original **file names** (e.g. `phenoDist_rank_04142022.csv`), so
outputs can be diffed against the originals. The legacy container maps `Z:/` the same way via a path overlay. Note the
`LOG/logFile*.txt` side effect of `printLog`, which is kept and redirected.

Notable inputs referenced but **not available**: `gene_key.csv`, `genes_wb_info3b.csv`, `EMBD_gene_info.csv`,
`ortholist_results_embd.csv`, `manualScoringForMySQL_2020update.csv`, `manualScoringForSeabornB.csv`,
`EMBD500genepilot_high_prev.csv`, `IntegratedInteractome.csv`, `pilot_for_heatmapAllSens.csv`, `cfsp_for_heatmap.csv`,
the Yanai cluster files, and `fileLookup.SKIPPED`. Only `gene_key.csv` (derivable) and `SKIPPED` (G6) affect the
published-output path.

## 8. Performance notes (for Phase 3, nothing changed now)

- `getPADS`: 503² Python-level `getPAD` calls, each recomputing `getDist2Zero` twice and allocating 276-length arrays,
  across 4 parameter sets. This is vectorizable to one masked-matrix operation.
- `getCSIS`/`getNBelow`: for each (i, j), a linear scan of `distAll` plus `r2 in rcounted` list membership. That is
  roughly O(n⁴) ≈ 6×10¹⁰ elementary operations for n=503, and is almost certainly the dominant cost. CSI is written to the CSV but is
  not used by the published rankings; it still has to be reproduced because it is in the CSV.
- `output_all_best_matches_*`: re-reads all four CSVs (≈250k rows each) **inside** the 503-iteration loop, and builds
  `RNAiClass` objects it never uses.
- `getGeneCommonName`: re-reads `gene_key.csv` on every call (≈10 per seed).

## 9. What I need from you before Phase 1

1. **Reference files** (A1): copy `all_best_matches_05192022.xlsx`, `EMBD0001_0504_…_052422.csv`,
   `README_Zenodo_EMBD_2023.md`, and `Index_EMBD_ID_Gene_names.xlsx` into `data/` (or tell me where they are). Also
   confirm that `data/` should be git-ignored or committed (licensing/size of the Zenodo files).
2. **Approve the layout** (§1.2): originals stay in place, unmodified; the new package goes under `src/embd/`.
3. **Approve the stub/shim approach** for the legacy container: original files untouched, with missing and
   out-of-scope modules stubbed on `PYTHONPATH`.
4. **Questions for the lab** (can run in parallel): A6 (why 10 rows are full precision; were values float32 at run time?),
   B10/P5 (which conditions are in the published universe), R2 (provenance of `controlOrigin`/`rnaiStd`), G6 (contents
   of `SKIPPED`), step-3 run order (§2), and the umap-learn version used in 2022 (C1).

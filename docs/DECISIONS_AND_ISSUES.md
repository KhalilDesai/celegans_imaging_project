# Decisions and Issues Log

Running list of factual findings, theoretical/scientific concerns, and design decisions
encountered while porting the EMBD downstream analysis. Policy for the first pass
(per maintainer, 2026-10-08): **port as faithfully as possible, bug-for-bug, making minimal
design decisions even where the original code looks wrong.** Anything that looks wrong is
recorded here instead of fixed. Fixes happen only after explicit approval and in their own
commits.

Status tags: `OPEN` (needs decision/info), `FAITHFUL` (kept as-is on purpose),
`DECIDED` (decision made, with who/when), `VERIFY` (claim to be confirmed empirically in Phase 1).

---

## A. Data / inputs

| # | Status | Finding |
|---|--------|---------|
| A1 | OPEN | Only `All_genes_mySQL_table.xlsx` is available (on the Desktop, not in the repo's `data/`). The reference outputs `all_best_matches_05192022.xlsx`, `EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv`, `README_Zenodo_EMBD_2023.md`, and `Index_EMBD_ID_Gene_names.xlsx` are missing. They're needed for Phase 1. |
| A2 | FACT | The sheet has 514 rows: the control (`embd_id` 0, no `mysql_id`), EMBD 1–503 (complete), and 10 extra conditions: 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905. That is 513 genes + control, not 514 + control. |
| A3 | FACT | This copy of the sheet contains **no** `'NULL'` strings. Missing values are empty cells (6,716 NaN cells). The brief expected `'NULL'`. The file may have been re-saved by Excel (an Excel lock file `~$All_genes_mySQL_table.xlsx` was present). The loader will treat both `'NULL'` and empty as NaN. |
| A4 | FACT | No `version`, `prev_GLS`, or `prev_MS` columns, as expected. `RNAiClass.load_params_data` treats `version is None` as "go to pickle → reprocess embryos". The data layer must supply a non-None `version`. `prev_GLS`/`prev_MS` are only stored and logged, and never used by PAD. |
| A5 | FACT | The sheet has 158 `*_GLS`/`*_MS` columns. 154 of them map to `PARAM_NAMES` entries (`PARAM_NAMES` has 138 names, i.e. 276 possible columns). The other 4 (`EMB_*`, `abL1_*`) are metadata. `bY`/`bYErr` are absent for both strains. **All parameters used by every parameter set (optimal 17 GLS + 24 MS; all 58 + 50; GLS-only; MS-only) are present.** In the original, `db_utils_embryos.get_row` silently skips non-existent columns, and `split_params` fills them with NaN, so absent columns behave identically. |
| A6 | VERIFY | **Precision.** The `genes` table parameter columns were created as MySQL `FLOAT` (single precision; `db_utils_embryos.create_param_names`, line 90). 504 of 514 rows hold values with ≤6 significant digits (consistent with MySQL `FLOAT` text rendering). 10 rows hold full double precision: embd 0 (control), 32, 73, 109, 150, 159, 233, 295, 440, 493. The original code saw whatever MySQLdb returned at run time (April 2022). If the DB stored float32 and the dump shows 6-digit rendering, the dump may differ from what Python saw by up to ~1 float32 ulp. Expect PAD differences of ~1e-6 relative from this alone. Why those 10 rows are full precision is unknown (perhaps edited after a later recompute?), so ask the lab. |
| A7 | FACT | Some `embd_id`s have parameters stored as Python `int` in Excel (e.g. `0`, `125251000`). The loader must cast all parameter columns to float64. 3,385 exact zeros exist and are real values, not missing. |
| A8 | OPEN | The sheet's `top_hits_PAD` / `reciprocal_matches` columns come from a different run (per the brief) and must not be used as the reference. |
| A9 | OPEN | `gene_key.csv` (EMBD↔sequence name↔common name, read by `embdFunc.getGeneName/getGeneCommonName/getEMBDNumber`) is not in the data. The data layer will derive it from the sheet (`embd_id`, `ce_id`, `gene_name`) or from `Index_EMBD_ID_Gene_names.xlsx` once available. Original encoding is latin-1 (`makeUMAPs_V2.py:17`). |

## B. PAD / ranking algorithm (faithful-port notes)

| # | Status | Finding |
|---|--------|---------|
| B1 | FAITHFUL | PAD (`RNAiClass.getPAD`) = `max(0, 1 - d(a,b) / (d(a,0) + d(b,0)))`, with `d` from `embdFunc.getNDDistance`: Euclidean distance over coordinates that are non-NaN in both vectors, scaled by `sqrt(N / (N - i1 - i2))`. Here `N` = number of parameters in the active set, and `i1`/`i2` = coordinates present in only one vector. No overlap → `d = inf` → PAD 0. |
| B2 | FAITHFUL | Coordinates: `posND[i] = (param - controlOrigin[key]) / rnaiStd[key]`, with all weights = 1. The origin and norms are **hardcoded dicts** in `RNAiClass.setOrigin` (the origin is the control mean; the norms are the std across RNAi conditions). If a norm == 0, the coordinate is not divided. Distance-to-control treats the control as 0 on every *active* parameter. |
| B3 | FAITHFUL | "Optimal" parameter subset = `params_use_assign._group_rank_opt()` (dated 3/7/18, th=4). It **removes** parameters whose rank count is **> 4**, giving 17 GLS + 24 MS. The other published sets (GLS-only, MS-only, all) are produced by swapping the returned function inside `get_params_use()` (commented alternatives). The port will make this a parameter, and the default stays optimal. |
| B4 | FAITHFUL | Top hits (`PADassist.PADFilter`): `np.argsort(pads)[::-1]`, then skip index 0 (assumed to be self, PAD 1.0) and take the next `nTop`. Ties use numpy's default (unstable) quicksort/introsort order. If another condition ties with self at PAD 1.0, "skip self" could skip the wrong row. |
| B5 | VERIFY | **Tie order is library-dependent.** numpy ≥1.25 (x86 AVX-512) and 2.x (AVX2) use SIMD sort kernels for `argsort`, so tie order may differ from numpy 1.13. Phase 1 will measure whether any ties exist inside each gene's top 51. |
| B6 | VERIFY | `phenoDist_rank_*.csv` writes PAD via `'{}'.format(x)`. Under Python 2.7 this is `str()` = 12 significant digits; under Python 3 it is shortest-repr (up to 17). The published rankings were computed from the **12-digit-rounded** CSV values, so rounding can create or break ties. The faithful port must write `%.12g` to reproduce the CSV byte-for-byte. Confirm the numpy 1.13 `float64.__format__` behavior in the legacy container. |
| B7 | FAITHFUL | `phenoDistCalc.getPADS` sorts each gene's rows by distance (Python stable sort) and reverses them, which sets the row order in the CSV, and that order feeds argsort tie-breaking. |
| B8 | FAITHFUL | **PAD_RANK** is a derived MySQL cache: for EMBD 1–503, `PADassist.update_pad_rank` stores `PADFilter(i, data_opt, 51)` into columns `rank1..rank50` / `pad1..pad50`. The 51st hit is silently dropped because no column exists, and missing slots default to 0. `pad*` columns are `FLOAT`, so they are float32. It is read by `db_utils_genes.get_reciprocal_matches` and `update_top_hits_PAD`. The port regenerates it in memory from the optimal PAD CSV. To be faithful, it applies float32 rounding to the stored PADs before the `> 0.55` / `> 0.5` threshold comparisons. |
| B9 | FAITHFUL | Reciprocal match: seed's top 10 with stored PAD `> 0.55` (strict); keep a hit if the seed appears in that hit's top-10 ranks, **with no PAD threshold in the reverse direction**. `top_hits_PAD`: top 10 with PAD `> 0.5`. Behavior also depends on `movement_GLS/MS` and `WT_like_or_phenotype` (`update_column_reciprocal_matches`). |
| B10 | OPEN | Main loops cover EMBD 1–503 only (`range(1, 504)`). The 10 extra conditions (674…1905) appear only in hand-entered lists (`PADassist.custom_updates_notes`), so some separate run that isn't in the code included them. Reciprocal lookups for a hit outside 1–503 would raise `IndexError` (no PAD_RANK row). Decide whether the port's default universe is 1–503 (faithful) or all 513, and check against the reference files. |
| B11 | FAITHFUL | `RNAiClass.paramsUseGLS/MS` are built from dict iteration (order differs py2 vs py3). They are only used for membership tests and counting, never for vector order (`posND` is indexed by the `PARAM_NAMES` list). So dict-order changes do **not** affect PAD. |
| B12 | FACT | `self.pNamesUse = np.array([GLS_list, MS_list])` builds a ragged 2-element object array. numpy ≥1.24 raises `ValueError`. Port: `np.array([...], dtype=object)`, which behaves the same. |
| B13 | FAITHFUL | `getPAD`: when `d(a,0)+d(b,0) == 0`, the result is NaN and `pad > 0` is False, so it returns 0. This includes a gene compared with itself if its distance to control is 0. Kept. |

## C. Groups / UMAP / figures

| # | Status | Finding |
|---|--------|---------|
| C1 | FACT | UMAP was **not** run under Python 2.7. `groupAssesment.make_UMAP_content_files` (py2) saves `PAD_Matrix_msParams.npy` etc. (its docstring says "needs python 3.6 to make the UMAP itself"). `makeUMAPs_V2.py` (py3.6) loads `PAD_Matrix_042122_optParams.npy`, so the files were renamed by hand, and the `_msParams` suffix is stale. `umap-learn` is not in requirements.txt, and its version is unknown. |
| C2 | VERIFY | UMAP randomness: `np.random.seed(71076)` then `umap.UMAP()` with defaults and **no `random_state`**. With `random_state=None`, umap-learn draws from the global numpy RNG (so it is seeded) but runs parallel numba SGD, which is not bit-reproducible across runs or thread counts. Faithful port: same seed and defaults. Reproducible mode (proposal, needs approval): `random_state=71076`, which forces single-threaded mode and changes coordinates. Exact UMAP coordinates from 2022 are not recoverable without the original umap/numba versions. |
| C3 | FACT | UMAP input = 503×503 PAD matrix (rows as feature vectors, default Euclidean metric). Hardcoded group lists and colors are duplicated in `groupAssesment.make_UMAP_content_files` and `makeUMAPs_V2.py`. EMBD 77 and 45 each appear in two groups. |
| C4 | FAITHFUL | `get_all_PAD2EO` / `get_all_dist2EO` drop NaN/inf entries **then** reshape to `sqrt(size)²`. Harmless for PAD (never NaN/inf). For distances, any `inf` would make the reshape fail or silently misalign. |
| C5 | FACT | `fileLookup.SKIPPED` is imported in four places (`RNAiClass.py:1996`, `db_utils_genes.py:835`, `groupAssesment.py:1582`, `make_figures.py:2349`) but is not defined in `fileLookup.py`. Those functions raise `ImportError` as shipped. |
| C6 | FACT | Modules imported but absent from the repo: `dataValidation` (groupAssesment module level; PADassist module level), `SkeletonAnalysis` (Embryos), `GermlineClass` (make_figures, function level). Any module importing `dataValidation` at module level cannot be imported as shipped. |

## R. RNAiClass normalization constants / control

| # | Status | Finding |
|---|--------|---------|
| R1 | FAITHFUL | `setOrigin` uses origin = hardcoded `controlOrigin` and norms = hardcoded `rnaiStd`. `controlStd` is used only for `stdV` (not in PAD), and `rnaiMean` is unused. These dicts are **scientific constants** and will be copied verbatim (byte-checked) into the port. |
| R2 | OPEN | Provenance: `get_control_origin()` (per-embryo control mean/std from the `embryos` table) is *probably* the source of `controlOrigin`/`controlStd`, but this can't be verified without that table. No code in the repo produces `rnaiMean`/`rnaiStd`. `tailHead_MS`/`devHead_MS` origin values were hand-edited ("CONTROL ORIGIN HAS BEEN CHANGED … 04/12/22"). Ask the lab. |
| R3 | FACT | `write_controlAvg_mysql` writes the hardcoded `controlOrigin` into `genes` row id 2627 (EMBD0000). So the sheet's control row (embd 0, full precision) is likely a **copy of the hardcoded origin**, not an independent measurement. Phase 1 will check this. The port does not read the control row for PAD. |
| R4 | FACT | `ControlClass()` needs no embryo data as long as its row has a non-None `version`. `defineFolders`/`getEmbryos` are image-only and buggy (`os.listdir` order with `[:20]`, bare `except`). |
| R5 | FACT | `embdFunc.getNDDistance` uses a bare `raise` on a dimensionality mismatch. That is a `RuntimeError` in py3 (and `TypeError: exceptions must be old-style classes` in py2). Behavior is equivalent ("crash"). Port: `raise ValueError` with the same message. This is a message-only change, logged as DECIDED-minimal. |

## U. Utilities / packaging

| # | Status | Finding |
|---|--------|---------|
| U1 | DECIDED (proposal) | `myFunc` (cv2, tifffile, PIL at module level) and `myMath` (lmfit, myFunc) block an image-free import. In-scope needs are only `readCSV`, `getSig`, `getSubIndex`, `nan_helper`, plus `myFigure` (pure matplotlib). The port will copy these verbatim into image-free modules, and the image functions become `NotImplementedError("image analysis not yet ported")` stubs. |
| U2 | FACT | `cumulativeDistribs`, `paramDistribs`, and most of `make_figures` need per-embryo data (the `embryos` table, Embryo objects) or images, so they are out of scope. In-scope exceptions: `cumulativeDistribs.getAllPCC/load_pcc/get_pca/get_lda` (logic only; `get_pca` sets NaN→raw 0 **before** normalizing, giving −origin/std, kept FAITHFUL), `make_figures.show_500gene_heatmap`/`show_manual_heatmap_by_gene` (need CSVs not provided), and `plot_tIntVsPtArrest`/`plot_maxIntGLSVsPtArrest` (hardcoded data). |
| U3 | FACT | `make_figures` integer-division bugs (`rows = len/columns` → float in py3 → `plt.subplot` error) are all in image-dependent functions, so they are out of scope. |
| U4 | VERIFY | matplotlib removals affecting in-scope plots: `hist(normed=)`, `set_axis_bgcolor`, tick `'on'/'off'`, `tick.label`, `edgecolor=''`, `myFigure.save(transparant=)` (misspelled kwarg). Figures are compared visually or structurally, not byte-for-byte. Python-3 fixes to these APIs are rendering-only and are logged per call site in Phase 2. |
| U5 | FACT | Unit tests `UnitTests/nDimTests.py` are the only image-free tests. Several are broken by design (compare `nanmedian` with a `nanmean` implementation; `.all()` on a list), so they are not usable as a reference. |

## G. groupAssesment (closest group / auto_description)

| # | Status | Finding |
|---|--------|---------|
| G1 | FACT | `groupAssesment` cannot be imported as shipped: line 14 `import dataValidation` (missing module), and its only consumer (`optimizeByInteractome`) is not on the live path. Port: lazy import + `NotImplementedError`. |
| G2 | FAITHFUL | The live path: `get_all_closest_groups` (EMBD **37–503**, not 1–503) → `find_closest_groups` → `get_closest_group_description`. It writes `genes.auto_description` (MySQL) or a CSV `closest_groups{%Y%m%d}.csv`. **`makeDataframe=True` returns before the SQL block**, so `to_csv=True, to_SQL=True` never writes SQL. |
| G3 | FAITHFUL | Centroids = `np.nanmean` of `posND` over the 17 hardcoded manual groups (`get_manual_groups`). Members 77, 31, 264, 45, 235, 261, 422 belong to more than one group. |
| G4 | FAITHFUL | **Gene→centroid PAD differs from `RNAiClass.getPAD`:** `inf` distance is replaced by 0 (so no overlap gives PAD = 1, the opposite of getPAD's 0), and `PAD ≤ 0` becomes NaN (getPAD gives 0). It uses the global `get_params_use()` lists. Faithful. Flag to the lab as possibly unintended. |
| G5 | FAITHFUL | Assignment thresholds: keep groups with `PAD > max − 0.05·max` (strict). If top PAD `< 0.5`, it is "poor match" unless the group is `dev delay` (at any PAD, due to `or` precedence) or `wt` with PAD > 0.35. If top PAD `< 0.55` and not moving, only the top group is kept. Descriptions are joined with `" \| "`. Sorting uses `df.sort_values` (quicksort, unstable), so ties could reorder. |
| G6 | FACT | `get_all_closest_groups` checks `'EMBD{0}'.format(i) in SKIPPED` (no zero padding), whereas `db_utils_genes` uses `{0:04}`. `SKIPPED` is undefined anyway (C5). Need SKIPPED's content from the lab, or treat it as empty (DECISION needed). |
| G7 | FACT | Broken under py2 too (not on live path, ported as-is or stubbed): `getRNAiGroupsPheno` returns 11 values but callers unpack 8; `getPAD2EO` returns string arrays used numerically; `optimizeByInteractome` calls `optimizeByGain` with 1 arg; `getDist2EO` `isinstance(x, object)` is always True; `make_tSNE_all_data` is defined twice (identical); `test_two_sample_ks_test` skips column index 1 (`maxG_GLS`). |
| G8 | VERIFY | Randomness: `getRNAiGroupsPheno` seeds `np.random.seed(1)` and then calls `train_test_split` without `random_state` (global RNG). tSNE uses `TSNE(random_state=71076, method='exact')` but sklearn defaults changed: `init` random→pca and `learning_rate` 200→auto in 1.2; `early_exaggeration` 4→12 in 0.19 (0.18.1 pinned, so 4.0). Faithful port must pass the 0.18.1 defaults explicitly. Bit-exact tSNE is unlikely. `LogisticRegression(random_state=50)`: solver default liblinear→lbfgs in 0.22. `ks_2samp` default changed to exact for small n in scipy 1.3. |

## P. Paper (Green et al. 2024, Cell, PMC11166207) vs. code

Source: full text via PubMed Central (DOI 10.1016/j.cell.2024.04.012). The PMC text extraction drops
equations and some inline numbers (e.g. the reported PAD confidence threshold), so equation-level
checks need the PDF or the STAR Methods equations.

| # | Status | Finding |
|---|--------|---------|
| P1 | CONSISTENT | Paper: 58 Germ Layer (GLS) + 50 Morphogenesis (MS) parameters, and an optimized 41-parameter subset (17 GLS + 24 MS) found by greedy downhill climbing on 17 manually curated groups (91 knockdowns). Code: `_all()` = 58 + 50 and `_group_rank_opt()` = 17 + 24. ✔ |
| P2 | CONSISTENT | Paper: "if two or more embryos failed to move, then only non-moving embryos were considered". Code: `RNAiClass.asignUse` uses `len(NoMove) > 1`. ✔ (This happens upstream of the parameter table, so it's out of scope here, but the table already reflects it.) |
| P3 | CONSISTENT | Paper: distance over shared parameters with a sqrt penalty only for dimensions present in exactly one vector. PAD = distance ÷ (sum of both distances to control). Code: `getNDDistance`, `getPAD` (as `1 − ratio`, clamped at 0). ✔ (The "1 −" and the clamp should be confirmed against Eq. in the PDF.) |
| P4 | DIFFERS | **Reciprocal matches.** Paper: "both genes are among each other's top 10 PAD-ranked lists". Code (`get_reciprocal_matches`) adds a forward-direction threshold: only seed hits with stored PAD **> 0.55** are considered (float32-stored), and there is no threshold in reverse. The published lists come from the code, so the port keeps the code behavior (FAITHFUL). |
| P5 | DIFFERS | Paper: "500 genes". Code loops over EMBD 1–503, and the sheet has 513 genes. Presumably 3 of 1–503 are controls, duplicates, or excluded, and the 10 high-numbered IDs are the follow-up conditions (Integrator/FACT/energy pathways, "not in the 500 gene list"). OPEN: confirm with the lab which 500. |
| P6 | DIFFERS | Paper: UMAPs "based on the PAD values" (optimal 41, all GLS, all MS). Code: UMAP runs on the 503×503 PAD **matrix rows as feature vectors** with Euclidean metric, not with PAD as a precomputed distance. Faithful port keeps the code behavior. The paper's UMAP overlay shows 11 groups (code's `groupLists`), while optimization used 17 groups. |
| P7 | NOTE | The paper's four broad groups (Cell Fate, Morphogenesis, GLDI, Morphogenesis-GL) use ">2 SD from control mean" on per-embryo control distributions and KS tests. These likely need **per-embryo control parameters** (`embryos` table via `db_utils_embryos.get_all_params`), which are not in the gene-level sheet. If so, regenerating the broad groups is blocked on embryo-level data, not only on image analysis. To be confirmed in the module audit. |
| P8 | NOTE | Paper: optimization used a 70/30-style train/test split of the 17 groups ("at least three members per group"). Matches the `_group_rank_opt` docstring. The optimization procedure itself is not needed to reproduce published outputs (its result is hardcoded), and porting it is optional. |

## D. Environment / reproducibility

| # | Status | Finding |
|---|--------|---------|
| D1 | FACT | requirements.txt is a Windows dev-machine freeze (py2exe, pywin32, vboxapi, Django, PyQt4…). It pins `matplotlib==1.4.0rc2` (a release candidate) and omits `umap-learn`, `pyensembl`, and `lmfit`'s numba/asteval transitive pins. |
| D2 | FACT | The `genes`-table pickles fallback (`RNAi_Data/EMBD####`) and `check_version=True` call git (`GitPython`) on the *analysis repo* to stamp versions. The faithful port keeps `check_version=False` (the default) and supplies a fixed version string from the data layer. |
| D3 | FAITHFUL | `phenoDistCalc` uses `multiprocessing.Pool(6)` for PAD and CSI. Results are order-preserving (`pool.map`), so parallelism doesn't affect output. |

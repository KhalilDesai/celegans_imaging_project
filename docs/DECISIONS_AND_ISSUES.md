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
| A1 | RESOLVED | All reference files were pulled from Zenodo record 7761843 (CC-BY-4.0). They sit in member `Zipped_EMBD_files/EMBD_summary_readme_files.zip` of part `-010.zip`, fetched via HTTP range requests (~1.7 MB rather than the 19.4 GB record). They are committed in `data/` with provenance and checksums in `data/README.md`. |
| A2 | FACT | The sheet has 514 rows: the control (`embd_id` 0, no `mysql_id`), EMBD 1–503 (complete), and 10 extra conditions: 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905. That is 513 genes + control, not 514 + control. |
| A3 | FACT (corrected) | The sheet **does** contain literal `'NULL'` strings (as the brief said): 5,969 cells, of which 2,050 are in parameter columns and 3,919 in metadata columns. Another 747 cells are truly empty (only 4 of them in parameter columns). My first check missed the `'NULL'`s because pandas `read_excel` silently converts `'NULL'` to NaN (default `na_values`). openpyxl shows them. The Desktop copy is byte-identical to the Zenodo copy (SHA-256 `16112262…`). The converter (`embd_data.sheet`) maps both `'NULL'` and empty to SQL NULL, i.e. NaN in the parameter dicts, matching `db_utils_embryos.get_row` (None → NaN). One gene has `auto_description` = NULL, which matters for `update_WT_like_or_phenotype` (substring test on None). |
| A4 | FACT | No `version`, `prev_GLS`, or `prev_MS` columns, as expected. `RNAiClass.load_params_data` treats `version is None` as "go to pickle → reprocess embryos". The data layer must supply a non-None `version`. `prev_GLS`/`prev_MS` are only stored and logged, and never used by PAD. |
| A5 | FACT | The sheet has 158 `*_GLS`/`*_MS` columns. 154 of them map to `PARAM_NAMES` entries (`PARAM_NAMES` has 138 names, i.e. 276 possible columns). The other 4 (`EMB_*`, `abL1_*`) are metadata. `bY`/`bYErr` are absent for both strains. **All parameters used by every parameter set (optimal 17 GLS + 24 MS; all 58 + 50; GLS-only; MS-only) are present.** In the original, `db_utils_embryos.get_row` silently skips non-existent columns, and `split_params` fills them with NaN, so absent columns behave identically. |
| A6 | VERIFY | **Precision.** The `genes` table parameter columns were created as MySQL `FLOAT` (single precision; `db_utils_embryos.create_param_names`, line 90). 504 of 514 rows hold values with ≤6 significant digits (consistent with MySQL `FLOAT` text rendering). 10 rows hold full double precision: embd 0 (control), 32, 73, 109, 150, 159, 233, 295, 440, 493. The original code saw whatever MySQLdb returned at run time (April 2022). If the DB stored float32 and the dump shows 6-digit rendering, the dump may differ from what Python saw by up to ~1 float32 ulp. Expect PAD differences of ~1e-6 relative from this alone. Why those 10 rows are full precision is unknown (perhaps edited after a later recompute?), so ask the lab. |
| A7 | FACT | Some `embd_id`s have parameters stored as Python `int` in Excel (e.g. `0`, `125251000`). The loader must cast all parameter columns to float64. 3,385 exact zeros exist and are real values, not missing. |
| A8 | OPEN | The sheet's `top_hits_PAD` / `reciprocal_matches` columns come from a different run (per the brief) and must not be used as the reference. |
| A9 | OPEN | `gene_key.csv` (EMBD↔sequence name↔common name, read by `embdFunc.getGeneName/getGeneCommonName/getEMBDNumber`) is not in the data. The data layer will derive it from the sheet (`embd_id`, `ce_id`, `gene_name`) or from `Index_EMBD_ID_Gene_names.xlsx` once available. Original encoding is latin-1 (`makeUMAPs_V2.py:17`). |
| A10 | FACT | Reference universe = **EMBD 1–503**. `all_best_matches_05192022.xlsx` has 503 seeds × top 10, no self-hits, and no hit IDs > 503. The 052422 CSV covers embd 1–503, and none of its hit names are the 10 extra conditions. So the faithful pipeline universe is `range(1, 504)`, exactly as the code loops. |
| A11 | FACT | `all_best_matches_05192022.xlsx` has **three** parameter sets (optimal, GLS, MS: 12 columns). That matches the 12-field `fd.write` in `output_all_best_matches_to_csv_combined`; the all-params columns are built in memory but never written. PADs carry **12 significant digits** (e.g. `0.811731716522`), i.e. py2 `str(float)`, which confirms B6. There are no duplicate PAD values within any seed's top 10. |
| A12 | FACT | The 052422 CSV went through Excel: gene `apr-1` appears as `1-Apr` in some cells (date auto-conversion), and the file is latin-1. Hit names `C01G5.2`, `cco-2`, `pacs-1` are not in the 2023 index, so gene names changed between the 2022 `gene_key.csv` and the 2023 index. Comparisons will map names→EMBD via the index plus an explicit alias table. |
| A13 | FACT | The sheet's `auto_description` equals the 052422 CSV for all 503 genes. `top_hits_PAD` differs in 246/503 rows and `reciprocal_matches` in 92/503 (cf. A8). |
| A15 | FACT | The sheet also lacks `WT_like_or_phenotype`, which `update_WT_like_or_phenotype` writes and `update_column_reciprocal_matches` reads (found when the legacy `derived` step failed with "no such column"). The converter adds it as NULL. The faithful pipeline computes it before reading it (L7). |
| A16 | FACT | The 2023 index also has an Excel date-mangled gene name: EMBD 87 (K04G2.8) `Gene Name` is the date 2022-04-01 (= `apr-1`). The gene-key reconstruction converts date cells back to `mon-day` (`apr-1`). Found because the first legacy `derived` run printed `2022-04-01 00:00:00` as a hit name. |
| A17 | OPEN | The paper's 500 genes are EMBD 1–500 (maintainer, Q2), but the code and every published output use EMBD 1–503. EMBD 501 (mex-6), 502 (lin-66) and 503 (cuti-1) are seeds and hits, and 501 and 503 belong to manual groups that set the group centroids, auto-descriptions and UMAP overlays. The faithful port keeps 1–503 because removing them would change published results. Restricting to 1–500 would be a deliberate, non-faithful option (a pipeline parameter) needing maintainer approval. |
| A14 | FACT | The two index files (`Index_…xlsx`, `A Index_…xlsx`) have identical cell contents (514 rows; columns `EMBD ID`, `C. elegans ID`, `Gene Name`) and differ only in xlsx metadata. EMBD 1903 is `hmg-3/hmg-4` and shares C. elegans ID `C32F10.5` with EMBD 3 (`hmg-3`). |

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
| B10 | RESOLVED (A10) | Main loops cover EMBD 1–503 only (`range(1, 504)`). The 10 extra conditions (674…1905) appear only in hand-entered lists (`PADassist.custom_updates_notes`), so some separate run that isn't in the code included them. Reciprocal lookups for a hit outside 1–503 would raise `IndexError` (no PAD_RANK row). Decide whether the port's default universe is 1–503 (faithful) or all 513, and check against the reference files. |
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

## L. Legacy baseline (Phase 1): how the original code is run

| # | Status | Item |
|---|--------|------|
| L1 | DECIDED | The original Python 2.7 files are run **unmodified**. The repo is mounted read-only into `embd-legacy:py27` (python:2.7.18-buster, amd64, pinned by digest). Missing or out-of-scope imports are satisfied by stubs ahead of the source on `sys.path` (`legacy_baseline/shims`). |
| L2 | DECIDED | **matplotlib 1.4.3** instead of the pinned `1.4.0rc2`: rc2 is not on PyPI, and 1.4.0 fails to build against Debian buster's FreeType (setupext `TypeError`). All other in-scope pins match requirements.txt exactly (numpy 1.13.3, scipy 1.0.0, pandas 0.21.1, scikit-learn 0.18.1, seaborn 0.8.1, networkx 1.11, lmfit 0.9.1, GitPython 2.1.8). matplotlib only affects figures, never PAD or rankings. |
| L3 | DECIDED | MySQL is replaced by a SQLite emulation (`src/embd_data/mysql_emulation.py`). Inferred MySQL semantics: double-quoted string literals; `SHOW COLUMNS`; FLOAT columns read back as `'%.6g' % float32(x)` (inferred from the dump, A6). Defaults: `genes` served as-is from the sheet, `pad_rank`/`pair_lookup` emulate FLOAT. Column types of the `genes` table are inferred (all-numeric → FLOAT, `id`/`embd_id` INTEGER, else TEXT). |
| L4 | DECIDED | Names absent from the public repo are supplied at run time: `varLookup.MYSQL_DB_CREDENTIALS=None`, `db_utils_embryos.mySQLpath/login` (dummies), `fileLookup.SKIPPED=[]` (Q4), `db_utils_genes.conn` (the `__main__`-level global its `update_*` functions use). |
| L5 | DECIDED | The parameter set is selected by replacing `params_use_assign.get_params_use` with `_group_rank_opt`/`_gls_params`/`_ms_params`/`_all` **before** `RNAiClass` is imported, one process per set. This is equivalent to the authors editing the `return` line. |
| L6 | DECIDED | `update_top_hits_PAD` (hand-edited list `[1406, 674, …]`) and `update_column_reciprocal_matches` (`range(310, 311)`) are run as verbatim copies with the loop set to EMBD 1–503. Python-2 `print(...).format()` is written as `print(....format())` (output text identical). `update_WT_like_or_phenotype` and `get_all_closest_groups` are called as-is. |
| L7 | DECIDED | Run order of the derived columns (Q5 open): `auto_description` → `WT_like_or_phenotype` → `top_hits_PAD` → `reciprocal_matches`. |
| L8 | DECIDED | Performance-only: in `bestmatches`, `PADassist.readCSV` is memoised (the original re-reads four ~253k-row CSVs on every one of 503 iterations). The function is pure, and its result is only read, so outputs cannot change. Disable with `EMBD_LEGACY_CACHE_READCSV=0`. |
| L9 | DECIDED | Date-stamped outputs (`Z:\\phenoDist_rank_<today>.csv`, `all_best_matches_<today>.csv`, `hits_matrix_v2_<today>.csv`, `*_msParams.npy`) are renamed to the names the downstream readers expect (`phenoDist_rank_04142022.csv`, `…GLS_ONLY05102022.csv`, …) or to `*_legacy.*`. |
| L10 | DECIDED | `gene_key.csv` is reconstructed from the 2023 index (A9, A12). Gene-name strings in legacy outputs can therefore differ from the 2022 originals for renamed genes. IDs are unaffected. |
| L11 | FACT | The legacy run is amd64 under emulation (Docker Desktop on Apple silicon), so floating-point semantics match the original x86-64 Windows workstation more closely than native arm64 would. |

## X. Phase 1 experiment: legacy code vs. published `all_best_matches_05192022.xlsx`

Run 2026-10-08 with the unmodified Py2.7 code (L1–L11) on the Zenodo sheet. Full numbers: `docs/reports/phase1_comparison.json`.

| Set | Exact ordered top-10 | Mean top-10 overlap | max abs dPAD | median abs dPAD |
|---|---|---|---|---|
| optimal (41) | 486/503 | 0.9984 | 0.0665 | 0 |
| GLS (58) | **503/503** | 1.0000 | 1.1e-16 | 0 |
| MS (50) | 487/503 | 0.9986 | 0.0228 | 0 |

| # | Status | Finding |
|---|--------|---------|
| X1 | FACT | **GLS reproduces bit-for-bit** (all 5,030 published PADs exactly, all orderings). This shows that the pipeline, the parameter selection, the hardcoded normalization, the 12-digit CSV text, and PADFilter ranking are reproduced faithfully. |
| X2 | FACT | The optimal-set differences > 1e-6 involve mainly **EMBD 1 (teg-4) and 478** (42 of 72 pairs). The rest are ~1e-6. A further 1,638 pairs differ by 1e-9–1e-6, spread over many genes (float32-scale changes). The MS-only differences involve **EMBD 181–185, 191, 363** (MS parameters only). |
| X3 | INFERENCE | The sheet (dumped 2023-03-21) is a **slightly different snapshot** from what each reference run read. The optimal file is dated 04/14/2022 and the GLS/MS files 05/10/2022. `fileLookup`'s note mentions "fixed arx-1 and arx-3 and fixed headInt params as of 4/14/22", and `RNAiClass` notes "CONTROL ORIGIN HAS BEEN CHANGED FOR tailHead and devHead as of 04/12/22". The genes that differ were most likely reprocessed between those dates and the dump. **The code is not the cause** (X1). → New lab question Q10. |
| X4 | DECISION NEEDED | Golden reference for regression tests: (a) the legacy outputs computed from the sheet (exactly reproducible), with the published file as a secondary check that allows the documented differences in X2; or (b) wait for the lab to provide the 2022 snapshot. Recommendation: (a). |
| X5 | FACT | Performance of the legacy run on 4 emulated cores: phenoDist ≈ 44–48 min per parameter set; `output_all_best_matches_to_csv_combined` **6.7 h** even with the readCSV memoization (L8). The cost is in `PADFilter` rebuilding `np.array(data)` (253k×7 mixed → strings) on each of ~4,000 calls. |
| X6 | FACT | **Published `auto_description` was produced by a newer `groupAssesment.py` than the public one.** 59 published descriptions contain "Embryos exhibit a variable point of arrest in the morphogenesis reporter strain. Embryos are mostly normal in germ layer strain.", and that text exists nowhere in the public code. The public `get_group_description` has only the germ-layer counterpart (`gMix`). The published run therefore had at least one extra manual group (and possibly different group lists). The legacy public code gives a different description for 50 of 467 recomputed genes (EMBD 37–503). The published column is **not reproducible** from the public code, so the faithful port targets the public code's output. `WT_like_or_phenotype` and `reciprocal_matches` depend on `auto_description` (the "slow pace" → Dev. Delay rule), so those differences propagate. → Q11. |
| X7 | FACT | Derived columns vs the published 052422 CSV (legacy public code, after A15/A16 fixes): `top_hits_PAD` identical for 488/503 genes, and **all 15 differences involve EMBD 1/478** (X2). `reciprocal_matches` identical for 475/503: 10 differences involve EMBD 1/478, 6 follow from `auto_description` differences (X6), and 12 are classification differences. In those 12, the published file says "no arrest phenotype" but the public code lists hits: 9 genes have `movement` 0 in one strain (public rule: either strain 0 → Phenotype), and 3 are "slow pace" (public rule → Dev. Delay). The published run used a different WT-like rule or run order (Q5, Q11). Name aliases learned between the 2022 and 2023 names (unambiguous): `cco-2`=`cox-5A`, `C01G5.2`=`prg-1`, `pacs-1`=`paic-1`. |
| X8 | FACT | Legacy wall-clock timings are **not** a valid performance baseline: the identical `bestmatches_combined` step took 6.7 h on the first run and 40 h on the rerun (laptop sleep/contention under amd64 emulation). Phase 3 will re-measure under controlled conditions (no sleep, a fixed workload subset, repeated runs). |

## P2. Python 3 port (Phase 2): decisions

| # | Status | Item |
|---|--------|------|
| P2-1 | DECIDED | Module and function names mirror the originals (`src/embd/RNAiClass.py`, `PADassist.py`, …), so each ported function can be diffed against its source. `groupAssesment.py` and `_rnai_constants.py` are **generated** from the original files by scripts in `tools/`, and the edits are listed in those scripts. |
| P2-2 | DECIDED | The parameter set is selected with `EMBD_PARAMSET` (opt/gls/ms/all), read by `get_params_use()`. The environment variable is inherited by multiprocessing workers (macOS spawns them), and each set runs in its own process as in the original (L5). |
| P2-3 | DECIDED | Python-2 number rendering in written files is reproduced by `_compat.py2_str` (`'%.12g'` plus `'.0'` for integral values), as verified against every sampled field of the golden CSV. Used in the phenoDist CSV writer, the `insert_row/update_row` SQL values, and the `str(val)[0:5]` truncation in hits_matrix_v2. |
| P2-4 | DECIDED | **Output-neutral speed-ups, per maintainer request (2026-10-16):** `PADassist.readCSV` is memoised on (path, mtime, size). `PADFilter`/`CSIFilter`/`distFilter` find a seed's rows through a label→row-indices index built once per list, instead of `np.where(np.array(data)[:, 0] == name)` on every call; `tests/test_port_units.py` tests the equivalence. A first version cached the whole `np.array(data)`, which is ~230 MB per file in numpy 2 because the mixed rows become fixed-width unicode, and the best-matches step was **OOM-killed** on this 8 GB machine (4 GB held by the Docker VM). It was replaced by the label index. |
| P2-5 | DECIDED | Already-broken original functions on no output path (`phenoDistCalc.optimizePAD`, `getCSIPAD`; most of `groupAssesment`'s experiment functions) are ported as-is or listed as not yet ported. They are never silently "fixed". |
| P2-6 | DECIDED | `PADassist.get_data_from_mysql` uses SQLite string syntax instead of MySQL `CONCAT/LPAD` (the emulation does not translate functions). It is not used by the pipeline, which reads the CSV (`get_data`). |
| P2-7 | FACT | Port timing (Apple-silicon laptop, native): phenoDist ≈ 56–86 s per parameter set, against 44–48 min for the legacy run under amd64 emulation. |
| P2-8 | FACT | macOS set the `hidden` file flag on the venv's `.pth` files after creation (the project sits in `~/Desktop`, possibly iCloud-synced), and CPython ≥ 3.12.8 skips hidden `.pth` files, so `import embd` failed mid-run. Fixed with `chflags nohidden`; `pipeline.py` now sets `PYTHONPATH` for its subprocesses. The production container (Phase 4) is unaffected. |

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

## Q. Open questions for the Oegema/Green lab

Collected 2026-10-08. Each links to the finding that raised it. Answers will be recorded here and the
corresponding item updated.

| # | Question | Raised by | Answer |
|---|----------|-----------|--------|
| Q1 | Why do 10 rows of the `genes` dump (embd 0, 32, 73, 109, 150, 159, 233, 295, 440, 493) hold full double-precision parameter values when all others have ≤6 significant digits? Were the parameter columns `FLOAT` (float32) when the April/May 2022 PAD files were generated? | A6 | **Answered 2026-10-26 (maintainer):** not a problem. Truncating everything to 6 significant figures would almost certainly be fine, but it is **not applied**: it would change PADs for those 10 genes against the frozen baseline for no benefit. It remains an available option. |
| Q2 | Which conditions make up the published universe? The paper says 500 genes, the code loops over EMBD 1–503, and the sheet has 513 genes (1–503 plus 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905). How were the 10 extra conditions added to the PAD files and best-match lists? | B10, P5 | **Answered 2026-10-26 (maintainer):** about 1,900 genes will eventually form the full set. EMBD **1–500** are the 500 of the first-pass paper. **But** the published outputs and the code also include **EMBD 501–503** (mex-6, lin-66, cuti-1): they are seeds and appear as hits in 51 of the 5,000 rows for genes 1–500, and two hardcoded manual groups contain them (501 in lMEX, 503 in gRupt). The faithful pipeline keeps 1–503; see A17. |
| Q3 | Where do the hardcoded `controlOrigin`, `controlStd`, `rnaiMean`, and `rnaiStd` dicts in `RNAiClass.setOrigin` come from (date, code, embryo set)? Was `get_control_origin()` the source of `controlOrigin`? What changed for `tailHead_MS`/`devHead_MS` on 04/12/22? | R1, R2 | — |
| Q4 | What did `fileLookup.SKIPPED` contain? It is imported in four places but absent from the repo. | C5, G6 | — |
| Q5 | In what order were the `genes`-table derived columns populated (`update_WT_like_or_phenotype`, `update_top_hits_PAD`, `update_column_reciprocal_matches`, `update_automated_description`)? `WT_like_or_phenotype` reads `auto_description`, so the order matters. | PORTING_PLAN §2 | — |
| Q6 | Which umap-learn (and numba, numpy) versions produced the 2022 UMAPs, and with how many threads? | C1, C2 | — |
| Q7 | Can you share `gene_key.csv` and `dataValidation.py` / `SkeletonAnalysis.py` (for completeness of the archive)? | A9, C6 | — |
| Q8 | Is the extra forward-direction PAD > 0.55 threshold in reciprocal matching intended (the paper describes only mutual top-10 membership)? | P4, B9 | — |
| Q9 | Is the gene→group-centroid PAD quirk intended (no parameter overlap gives PAD = 1 because inf distance is replaced by 0)? | G4 | — |
| Q10 | The published sheet differs from the inputs of the April/May 2022 PAD runs for EMBD 1, 478 (optimal set) and 181–185, 191, 363 (MS set), plus ~1e-7-level differences elsewhere. Were these genes reprocessed after 05/2022 (e.g. the head-parameter / arx fixes)? Is a dump of the `genes` table as of 05/2022 available? | X2, X3 | **Answered 2026-10-26 (maintainer):** not important; record as a discrepancy and use the latest data (the Zenodo 2023 sheet) as the baseline. This is already the case (X4 option a), and the differences stay documented in X2/X3. |
| Q11 | The published `auto_description` (052422) uses a group description absent from the public `groupAssesment.py` ("variable point of arrest in the morphogenesis reporter strain …", 59 genes). Can you share the version of `groupAssesment.py` (`get_manual_groups`, `get_group_description`) used on 05/24/2022? | X6 | — |

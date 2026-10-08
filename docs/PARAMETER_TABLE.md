# Per-gene parameter table: the image-analysis ↔ downstream interface

The downstream package (PAD similarity, rankings, matches, groups, UMAP, figures) takes exactly one
input: a **per-condition parameter table** with the columns of the lab's MySQL `genes` table. The
published instance is `data/All_genes_mySQL_table.xlsx` (Zenodo 7761843). A future port of the
image-analysis code (`embryos_analyze/Embryos.py` → `RNAiClass.calcAverageParams`) only has to produce
this table to plug in.

## 1. Rows

- One row per RNAi condition, keyed by **`embd_id`** (int). `embd_id = 0` is the control. In the 2022 data
  the conditions are EMBD 1–503 plus 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905. The published
  pairwise outputs use **1–503 only** (DECISIONS A10).
- **`id`** (int): the MySQL primary key. The original code looks rows up by `id` after resolving `embd_id`
  (`db_utils_genes.get_id_from_embd`). The dump calls it `mysql_id`. The control has none, and the data layer assigns
  2627 (DECISIONS R3).

## 2. Parameter columns

Name pattern: **`{p}_GLS`** (Germ Layer strain) and **`{p}_MS`** (Morphogenesis strain), for each `p` in
`varLookup.PARAM_NAMES` (138 names). Names that already carry a strain tag keep it, so the column becomes e.g.
`CoM0Rtd_GLS_GLS`, `MoI1GstdTail_MS_MS`.

**Value semantics.** Each value is the **median, over the condition's usable embryos of that strain**, of
the per-embryo parameter (`RNAiClass.calcAverageParams`):
- Usable embryos are the non-moving ones if more than one embryo failed to move, otherwise all non-dead embryos
  (`RNAiClass.asignUse`; paper STAR Methods).
- An embryo contributes only if `Embryo.checkGoodParam(p)` accepted its fit for that parameter.

**Required columns.** Only parameters in an active parameter set enter PAD. The union of the four published sets
is required (58 GLS + 50 MS names, defined in `development/params_use_assign.py`). Every other `PARAM_NAMES`
column is optional, and an absent column behaves exactly like an all-NaN column (DECISIONS A5).

| Family (`p`) | Meaning | Unit / scale |
|---|---|---|
| `aG, aR, aY` | sigmoid plateau of nuclear counts (GLS) or total intensity (MS), per channel G/R/Y | nuclei (GLS); summed camera counts above threshold (MS) |
| `bG, bR, bY` | sigmoid baseline | same as `a*` (0 for nuclear counts) |
| `mG, mR, mY` | sigmoid midpoint | time points (20 min) after the embryo's `t0`, time-scaled (from code: `Embryo.dt = 1`) |
| `sG, sR, sY` | sigmoid width | time points |
| `rG, rR, rY`, `*Err` | sigmoid shape term, fit standard errors | as the parent parameter |
| `Gtail, Rtail, Ytail` | mean of the last 5 (G/R) or 7 (Y) time points | nuclei |
| `residsEarly*`, `residsLate*` | mean absolute residuals before/after truncation (signed for late) | nuclei |
| `maxG, maxR` | max of the total-intensity curve after the first 10 time points | camera counts |
| `mRmG, mYmG, mYmR` | difference of sigmoid midpoints | time points |
| `fracR, fracG, fracY` | fraction of nuclei of each colour | dimensionless |
| `tMove` | movement onset, `(tMove − t0)·tScale` | time points |
| `movement` | 1 if sustained movement detected, else 0 (per embryo). The gene-level median is a fraction-like value. | dimensionless |
| `tScale` | time-scaling factor aligning the embryo to control | dimensionless (≈1) |
| `CoM{0,1}{G,R,Y}{td,scale,avgRes,stdTail}_{GLS,MS}` | centre-of-mass distance feature (axis 0/1, channel): `td` time of deviation from control (capped at the control movement time), `scale` y-scaling vs control, `avgRes` mean residual after `td`, `stdTail` std after `td` | `td`: time points; `scale`: dimensionless; `avgRes`/`stdTail`: feature units (distances normalised by embryo radius / half-length) |
| `MoI{0,1}{G,R,Y}{…}_{GLS,MS}` | normalised moment-of-inertia (radius of gyration) shape feature, same four sub-parameters | as above |
| `{a,m,s}SigHead`, `{a,m,s}GaussHead`, `maxHead, tailHead, scaleHead, devHead` | head-enclosure feature (MS end-on epidermal area) | area fraction / time points (inferred) |
| `scaleLength, devLength, tailLength` | neuronal-length feature (MS) | normalised length (inferred) |

Units marked *inferred* come from the paper's STAR Methods and the variable names, not from an
explicit unit in the code. A producer should still document its own units. The PAD depends only on
`(value − control origin) / RNAi std` (DECISIONS B2), so any change of units must be applied consistently to
the hardcoded normalisation constants in `RNAiClass.setOrigin`.

## 3. Missing-value convention

- **NaN (empty cell or the string `'NULL'`) means the parameter is unavailable for that condition**, e.g. a
  marker was not expressed, no embryo passed `checkGoodParam`, or the strain was not imaged. It is **not** 0.
- `0` is a legitimate measured value (3,385 zeros in the 2022 table; DECISIONS A7).
- PAD handles NaN through overlap and a penalty: coordinates present in only one of the two vectors inflate the
  distance by `sqrt(N / (N − i1 − i2))` (DECISIONS B1).

## 4. Precision

Store float64. The 2022 values passed through MySQL `FLOAT` (single precision) and were rendered with 6
significant digits (DECISIONS A6, Q1). A producer that writes float64 will give PADs that differ slightly
(~1e-6 relative) from the published ones.

## 5. Non-parameter columns read downstream

| Column | Used by | Notes |
|---|---|---|
| `movement_GLS`, `movement_MS` | reciprocal matches, WT-like call, closest group | compared to exactly 0 and 1 |
| `auto_description`, `WT_like_or_phenotype`, `top_hits_PAD`, `reciprocal_matches` | **written** by the pipeline (step 3), and `auto_description` / `WT_like_or_phenotype` are read back | inputs are overwritten. Their initial values matter only for EMBD 1–36 `auto_description` (not recomputed; DECISIONS G2) |
| `man_group`, `man_GS*`, `man_MS*` | classifier / KS experiments only | not on the published-output path |
| `gene_name`, `ce_id` | labels only (via `gene_key.csv`) | — |

## 6. Data-layer bookkeeping (not science)

`version` (non-null string), `prev_GLS`, `prev_MS` (may be NaN). `RNAiClass.load_params_data` reads them, and
a NULL `version` makes the original code fall back to pickles and embryo reprocessing. The data layer
always supplies `version = "zenodo-7761843"` (DECISIONS A4).

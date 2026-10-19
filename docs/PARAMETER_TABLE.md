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

## 3. Missing-value convention: 0 vs NULL

The convention is not documented by the authors, but the original code applies it consistently.
**0 means "measured, and the answer is zero/absent". NULL (NaN) means "undefined, or could not be
measured".** A producer of this table must reproduce it, because PAD treats the two very differently.

### 3.1 Per embryo (`embryos_analyze/Embryos.py`)

Every parameter starts as NaN (l. 119) and is overwritten only when it can be computed.

- **Marker never switches on → amplitudes 0, timings NaN.** Instead of fitting a sigmoid, the code assigns:

  | Strain, trigger | Set to 0 | Set to NaN | Other |
  |---|---|---|---|
  | GLS `setSigmParams` (l. 1705): ≤ 5 time points with > 20 nuclei in that colour | `a*`, `s*`, `*tail`, `*tailErr`, `aErr`, `sErr` | `m*`, `r*`, `residsEarly*`, `mErr`, `rErr` | `residsLate*` = mean of the control curve |
  | MS `setSigmParams` (l. 776): total intensity `maxG` ≤ 3e8 or `maxR` ≤ 1.5e8 | `a*`, `b*`, `s*` and their errors, `rErr` | `m*`, `r*`, `mErr` | `movement` = 0, `tMove` = NaN |

  Rationale: "zero nuclei / no signal" is a measurement; "the midpoint of a curve that never rose" is undefined.
- **Shape/position curves (CoM*, MoI*; `setCM_MoIParams`, l. 380):**
  - ≤ 2 valid points → **sentinel values, not NaN**: `td = 0`, `scale = 1`, `avgRes = −mean(control curve)`, `stdTail = 0`;
  - the curve deviates only at its last points → `avgRes = 0`, `stdTail = 0`.
- **Derived differences** (`mRmG`, `mYmG`, `mYmR`) are computed only if both midpoints pass `checkGoodParam`.
  Otherwise they stay NaN.
- **Head-enclosure parameters (MS)** when the sigmoid+bump fit is not attempted: `mSigHead`, `sSigHead`,
  `mGaussHead`, `sGaussHead` = NaN; `aGaussHead` = 0. Since 04/07/2022, `tailHead`, `scaleHead`, `devHead` are
  computed in both branches. Before that they were 0/NaN (comments at l. 904–908), which is likely one reason the
  published sheet differs from the 2022 runs for some genes (DECISIONS X3).

### 3.2 Per condition (`RNAiClass.calcAverageParams`, `Embryo.checkGoodParam` l. 510)

The gene value is the median over usable embryos whose value passes `checkGoodParam`:
- NaN never passes;
- a positive value with an error passes if error/value < 0.5;
- 0 (or a negative value) always passes, because the 50% test only applies to positive values.

So a gene-level value is:
- **0** when its embryos genuinely measured zero;
- **NULL** when no embryo produced a usable value: the parameter was undefined, every fit was too uncertain, or
  no embryos of that strain were usable.

NaN is written to MySQL as `NULL` (`db_utils_embryos.update_row`).

### 3.3 Consequence for PAD

- **0 is a real coordinate.** After normalisation `(0 − control mean) / RNAi std`, it lies far from the
  control origin and contributes strongly to distances ("no endoderm" is a large phenotype).
- **NULL drops the dimension.** If only one of two genes has it, the distance is inflated by
  `sqrt(N / (N − i1 − i2))` (DECISIONS B1). If neither has it, it is ignored.
- There are 3,385 zeros in the 2022 table (DECISIONS A7).

### 3.4 In the published sheet

Missing values appear as the text `'NULL'` (2,050 parameter cells) and as truly empty cells (4 parameter cells).
Nothing in the code distinguishes them. Both load as NaN, as in the original MySQL path (DECISIONS A3).

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

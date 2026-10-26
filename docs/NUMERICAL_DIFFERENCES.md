# Numerical differences: Python 3 port vs. legacy (Python 2.7) baseline

Reference: `tests/golden/` (legacy outputs, DECISIONS X4 option (a)). Each entry gives what differs,
why, by how much, and the test tolerance with its justification. Entries marked **PENDING** need
maintainer approval (rule: a library upgrade that changes numerical output must be shown and approved).

## N1. Euclidean norm summation order (BLAS) — APPROVED 2026-10-16 (maintainer)

- **What:** `np.linalg.norm` (used by `embdFunc.getNDDistance`) computes `sqrt(dot(x, x))` through
  BLAS. Legacy: numpy 1.13.3 with its bundled x86-64 OpenBLAS (under amd64 emulation). Port:
  numpy 2.5 on its default BLAS. Different BLAS kernels sum the squares in a different order.
- **Size (optimal parameter set, all 503×503 pairs):**
  - per-gene ND vectors `posND`: **bit-identical** (all four sets);
  - `dist2zero`: up to 7.1e-15 absolute (≈1 ulp; relative ≈2e-16);
  - `phenoDist_rank_04142022.csv`: 10 of 253,009 rows differ, all in the **PAD** column, all in the
    12th significant digit (max abs 1.0e-12, max rel 6.4e-12). `dist`, `dist2WT`, `dims`, `CSI` text
    identical. Row order identical.
  - rankings: top-10 and top-51 lists identical for **503/503** seeds. None of the 10 pairs is in any
    seed's top 51, so `all_best_matches`, PAD_RANK, top hits and reciprocal matches are unaffected.
- **Full run (all four parameter sets):** 44 printed values in total differ across the four files
  (opt 10 PAD; GLS 6 PAD + 4 `dist`; MS 14 PAD + 4 `dist`; all 4 PAD + 2 `dist`). **Every one is
  exactly 1 unit in the 12th significant digit.** `dist2WT`, `dims` and `CSI` are text-identical, row
  order is identical, and top-51 rankings are identical for 503/503 seeds in every set.
- **Extension (APPROVED 2026-10-26, maintainer):** the approved rule is applied to the
  `dist` column as well (same cause). Tests allow ≤ 2 units in the 12th significant digit for
  `dist` and `PAD`, and exact text for every other column.
- **Original proposal (approved):** exact equality for row order, rankings and every non-PAD column; PAD
  compared numerically with `atol = 2e-12` (two units in the 12th printed digit: the legacy value
  itself was rounded to 12 significant digits, and a 1-ulp change in the underlying double can flip
  that last digit). The full-precision PAD/distance matrices (`pad_matrix_*.npy`) are compared with
  `rtol = 1e-13` (≈ 450 ulp head-room over the observed 1-ulp differences in distances; PAD is a ratio
  of distances, so error stays at this order).
- **Alternative (not recommended):** reimplement OpenBLAS-x86 `ddot` summation order in Python
  so that the 12-digit text is byte-identical. Fragile (kernel- and CPU-dependent), slow, and
  it buys nothing scientifically.

## N2. PAD_RANK stored PADs: legacy SQLite decimal parsing (harness artifact)

- **What:** the legacy driver stored PADs in SQLite as 12-digit text. The SQLite in the Python 2.7
  image (Debian buster, 3.27) converts decimal text to double with a routine that is not correctly
  rounded, so 10 stored values are 1 ulp off (e.g. `0.525694835345` → `0.5256948353449999`). The
  port's SQLite rounds correctly.
- **Effect:** none on outputs. The code reads these columns through MySQL FLOAT rendering
  (`'%.6g'` of float32), where both values are identical, so the reciprocal-match thresholds
  (> 0.55, > 0.5) see the same numbers. All 503×50 ranks are identical.
- **Tolerance:** ranks exact; PADs equal after `mysql_float_render` and within 1 ulp raw.

## N3. UMAP embeddings: not reproducible in the original either

- **What:** `makeUMAPs_V2` calls `umap.UMAP()` with no `random_state` (it seeds only the global numpy
  RNG). umap-learn then uses parallel, non-deterministic SGD. Three reruns of the **unmodified
  legacy** script (Python 3.6, umap-learn 0.5.3) all differ from the golden run.
- **Size (optimal set, 503 genes, 15-nearest-neighbour overlap with the golden embedding):**
  legacy reruns vs each other and vs golden 0.781–0.83 (mean 0.80); port (umap-learn 0.5.12, 4 runs)
  vs golden 0.779–0.79; port vs port 0.90. Trustworthiness (k=15, in PAD space): golden 0.969, legacy
  reruns 0.969–0.970, port 0.970–0.971. Procrustes disparity vs golden: legacy 0.023–0.055, port
  0.016–0.054.
- **Tolerance:** kNN-15 overlap with golden ≥ 0.75 and trustworthiness ≥ 0.96, i.e. inside the
  original code's own run-to-run variation (min 0.781 / 0.969), with a small margin.
- **Decision (APPROVED 2026-10-26, maintainer):** the port seeds UMAP (`random_state=71076`, which
  umap-learn runs single-threaded), so repeated runs give identical embeddings. They are still not equal
  to any 2022 embedding; the structural tolerance above still applies against the golden run. The
  original unseeded call stays available as `makeUMAPs_V2.run(..., random_state=None)` ("faithful mode").

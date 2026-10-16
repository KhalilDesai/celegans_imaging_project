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
- **Proposed tolerance:** exact equality for row order, rankings and every non-PAD column; PAD
  compared numerically with `atol = 2e-12` (two units in the 12th printed digit: the legacy value
  itself was rounded to 12 significant digits, and a 1-ulp change in the underlying double can flip
  that last digit). The full-precision PAD/distance matrices (`pad_matrix_*.npy`) are compared with
  `rtol = 1e-13` (≈ 450 ulp head-room over the observed 1-ulp differences in distances; PAD is a ratio
  of distances, so error stays at this order).
- **Alternative (not recommended):** reimplement OpenBLAS-x86 `ddot` summation order in Python
  so that the 12-digit text is byte-identical. Fragile (kernel- and CPU-dependent), slow, and
  it buys nothing scientifically.

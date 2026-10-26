# Randomness and determinism

What can make two runs of the pipeline differ, and what the port does about it.

| Source | Where | Handling |
|---|---|---|
| UMAP stochastic optimisation | `makeUMAPs_V2.run` | Seeded: `np.random.seed(71076)` (as the original) **and** `umap.UMAP(random_state=71076)` (maintainer decision, NUMERICAL_DIFFERENCES N3). umap-learn then runs single-threaded, so repeated runs on one platform give identical embeddings (`tests/test_port_units.py::test_umap_seeded_is_deterministic`). `random_state=None` restores the original unseeded, parallel, non-reproducible call. |
| Thread count in BLAS / numba | numpy dot products (distances, PAD), UMAP (numba) | The images set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMBA_NUM_THREADS` to 1 (Dockerfile, apptainer.def). |
| BLAS implementation (summation order) | `sqrt(x.dot(x))` in the distance | Not random, but platform-dependent: Apple Accelerate (native macOS) and OpenBLAS (the image) round the last bit differently. Effects are ≤ 1 unit in the 12th printed digit and change no ranking (N1, "Platform"). UMAP coordinates are sensitive to those bits, so compare UMAPs by structure (kNN overlap), not by value. **For byte-identical outputs, run the image** (same OS, libraries and BLAS on every machine of the same architecture). |
| `multiprocessing.Pool` | original phenoDist/CSI | No longer used by the port (row-vectorised, PERFORMANCE #2/#4); the originals used order-preserving `pool.map`, so it never affected results (DECISIONS D3). |
| Dictionary / set ordering | Python 3 | Insertion-ordered dicts; no output depends on set iteration order (checked against the golden outputs). |
| Dates in file names | original `time.strftime` names | The pipeline passes fixed file names (those of the legacy baseline). |

`makeUMAPs_V2` is the only place that draws random numbers in the pipeline (`myFigure` has a demo function
using `np.random.normal` that the pipeline never calls).

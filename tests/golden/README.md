# Golden reference outputs

Frozen outputs of the **original Python 2.7 code**, run unmodified by `legacy_baseline/` on
`data/All_genes_mySQL_table.xlsx`. They are the regression reference for the Python 3 port (maintainer
decision 2026-10-16, DECISIONS X4 option (a)). The published Zenodo files remain a secondary check,
with known, documented differences (DECISIONS X1–X7).

- `*.xz`: xz-compressed, exact bytes. `CONTENT.sha256` holds the hashes of the decompressed content.
- `MANIFEST.sha256`: hashes of the stored files, verified by `tests/test_golden_integrity.py`.
- `provenance.json`: images, package versions, input checksum, harness commit, legacy timings.
- UMAP files come from the Python 3.6 stage (`umap-learn 0.5.3`, a guessed version; DECISIONS C1/C2) and are
  **not** bit-reproducible (parallel SGD). They are compared structurally, not exactly.

Regenerate (≈ 2–3 days of compute under amd64 emulation; see `legacy_baseline/README.md`), then:

    python tests/golden/build_golden.py "build/legacy_work/Z:" build/umap_work tests/golden

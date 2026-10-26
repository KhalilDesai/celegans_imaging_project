# Questions for the Green / Oegema lab about EMBD_public

Context: we are porting the downstream part of the EMBD analysis (per-gene parameters → PAD
similarity → best matches, reciprocal matches, auto-descriptions, UMAP) from Python 2.7 to Python 3,
starting from the Zenodo data (doi:10.5281/zenodo.7761843) and the public GitHub code. Running the
original code unmodified on the Zenodo `All_genes_mySQL_table.xlsx` reproduces
`all_best_matches_05192022.xlsx` exactly for the Germ Layer parameter set (503/503 genes). For the
optimal and Morphogenesis sets it matches 486/503 and 487/503 genes. A few questions would help us
pin down the remaining differences. Internal references (A6, X2, …) point to our decisions log.

## Answered by the maintainer (2026-10-26)

- **Q1 (snapshot):** not important. We use the latest data (the 2023 Zenodo table) as the baseline and
  keep the differences documented.
- **Q2 (precision):** not a problem. Values are used as published.
- **Q3 (which 500):** EMBD 1–500 are the 500 genes of the first-pass paper (about 1,900 in the full set).
  Follow-up for the lab: the code and published outputs also include EMBD 501–503 (mex-6, lin-66,
  cuti-1), including in two manual groups. Were these intended to be part of the analysis?

## Data snapshot

1. **Which genes changed after the May 2022 runs?** The best-match lists differ only where EMBD 1
   (teg-4) and EMBD 478 (optimal set) or EMBD 181–185, 191 and 363 (Morphogenesis set) are involved.
   Many other pairs also differ by ~1e-7, which looks like float rounding. Were these genes reprocessed
   between May 2022 and the March 2023 table dump, perhaps with the head-parameter or arx-1/arx-3 fixes
   mentioned in the code comments? Is a dump of the `genes` table as of May 2022 available? *(X2, X3)*
2. **Precision of the stored parameters.** Most rows of the dump show 6 significant digits, consistent
   with the MySQL `FLOAT` columns the code creates. Ten rows (embd 0, 32, 73, 109, 150, 159, 233, 295,
   440, 493) have full double precision. Why do those rows differ? *(A6)*
3. **Which conditions are "the 500 genes"?** The paper says 500, the code loops over EMBD 1–503, and the
   table has 513 (1–503 plus 674, 703, 866, 885, 952, 1360, 1406, 1903, 1904, 1905). Which are in the
   published set, and how were the 10 extra conditions added to the PAD files? *(B10, P5)*

## Code that is not in the public repository

4. **Newer `groupAssesment.py`.** 59 published auto-descriptions (052422 CSV) use the text
   "Embryos exhibit a variable point of arrest in the morphogenesis reporter strain. Embryos are mostly
   normal in germ layer strain." That text does not appear anywhere in the public code. Could you share
   the version of `groupAssesment.py` (`get_manual_groups`, `get_group_description`) used in May 2022?
   *(X6)*
5. **`fileLookup.SKIPPED`.** It is imported in four places but is not defined in the public
   `fileLookup.py`. What did it contain? *(C5, G6)*
6. **Missing files and modules:** `gene_key.csv` (we rebuilt it from the Zenodo index; a few gene names
   have changed since 2022), `dataValidation.py` and `SkeletonAnalysis.py`. Could these be added to the
   archive? *(A9, C6)*

## How the published outputs were produced

7. **Run order of the `genes`-table updates.** `WT_like_or_phenotype` reads `auto_description`,
   and `reciprocal_matches` reads `WT_like_or_phenotype`. In what order were `update_automated_description`,
   `update_WT_like_or_phenotype`, `update_top_hits_PAD` and `update_column_reciprocal_matches` run?
   For 12 genes the published file says "no arrest phenotype" where the public code's rules
   (movement 0 in one strain → phenotype; "slow pace" → developmental delay) would list hits. *(Q5, X7)*
8. **UMAP environment.** Which versions of umap-learn, numba and numpy produced the published UMAPs, and
   with how many threads? The script calls `umap.UMAP()` without a `random_state`, so the embeddings
   cannot be regenerated exactly. *(C1, C2)*

## Methods / intent

9. **Hardcoded normalisation constants.** In `RNAiClass.setOrigin`, where do `controlOrigin`,
   `controlStd`, `rnaiMean` and `rnaiStd` come from (which code, which embryos, which date)? Was
   `get_control_origin()` the source? What changed for `tailHead_MS`/`devHead_MS` on 04/12/22? *(R1, R2)*
10. **Reciprocal-match threshold.** The paper describes reciprocal matches as genes in each other's top 10.
    The code additionally requires PAD > 0.55 in the forward direction. Is that intended? *(P4, B9)*
11. **Group-centroid PAD.** In `find_closest_groups`, when a gene shares no parameters with a group
    centroid, the distance becomes 0 and the PAD becomes 1 (a perfect match), the opposite of the
    gene-to-gene PAD, which gives 0. Is that intended? *(G4)*

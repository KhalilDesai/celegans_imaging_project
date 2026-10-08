# Reference data

These files are a verbatim copy of the summary files in the Zenodo record for Green et al. 2024. They are not modified.

- **Record:** Green RA, Khaliullin RN, et al. *Automated Profiling of Gene Function During Embryonic Development* (data).
  Zenodo, [doi:10.5281/zenodo.7761843](https://doi.org/10.5281/zenodo.7761843), published 2024-03-22.
- **License:** [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/). Credit: Green, Khaliullin, Zhao, Ochoa, Hendel,
  Chow, Moon, Biggs, Desai & Oegema (Ludwig Institute / UC San Diego).
- **Source location:** record file `Zipped_EMBD_files-20231206T131336Z-010.zip` → member
  `Zipped_EMBD_files/EMBD_summary_readme_files.zip` (1,724,858 bytes). It was fetched on 2026-10-08 with HTTP range requests, and
  the outer zip's CRC was verified by `zipfile`. The other record files (≈19.4 GB of per-condition movies, UMAP images and
  hits matrices) are not included.
- **Paper:** Green et al. (2024) *Cell* 187, [doi:10.1016/j.cell.2024.04.012](https://doi.org/10.1016/j.cell.2024.04.012).

| File | Role in this project |
|---|---|
| `All_genes_mySQL_table.xlsx` | **Input**: dump of the MySQL `genes` table (513 conditions + control). The boundary interface is described in `docs/PARAMETER_TABLE.md`. |
| `all_best_matches_05192022.xlsx` | **Primary reference output**: top 10 PAD matches per seed (EMBD 1–503) for the optimal, GLS-only and MS-only parameter sets. |
| `EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv` | **Reference output**: `top_hits_PAD`, `reciprocal_matches`, `auto_description` per gene (latin-1, has been through Excel; see DECISIONS A12). |
| `Index_EMBD_ID_Gene_names.xlsx`, `A Index_EMBD_ID_Gene_names.xlsx` | EMBD ID ↔ C. elegans ID ↔ gene name (identical cell contents). |
| `README_Zenodo_EMBD_2023.md` | Original Zenodo README. The `.txt` twin in the archive is byte-identical and not copied. |

`All_genes_mySQL_table.xlsx` is byte-identical to the copy the maintainer provided separately.

## SHA-256
```
928b48f0e94d766c631446be8f17ed6371a9a1c8a4d550663a07067285e39d3c  A Index_EMBD_ID_Gene_names.xlsx
16112262f185a0b3b3593676990896a9e7402b7340d4789d86550ee80fefaea6  All_genes_mySQL_table.xlsx
350cecd9937a106952ee36168947242bfbe3c5115a1ed793b4988d917a46b01b  Index_EMBD_ID_Gene_names.xlsx
f6ed4c423ebaa199e2fca0b8da0191d5c27cb1cc0ca7c01ce931aed966f85729  all_best_matches_05192022.xlsx
347dabaf02fb45d16d4cafd6a3353d9be53317373633830dacf1a8417b3fe3f1  EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv
8aa868cef46534273e9a0d163bdeb7eba820504227fe8beaf21d8a6263d53666  README_Zenodo_EMBD_2023.md
```

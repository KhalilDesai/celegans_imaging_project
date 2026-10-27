# Project notes: index

Every note about the Python 3 port, in one list. Each file has one job; add new notes as their own file and
list them here.

## Start here

1. `HANDOFF.md`: current status, what is left, working agreements, how to set up a machine.
2. `../notebooks/walkthrough.ipynb`: one gene (EMBD0087, *apr-1*) through the whole pipeline, with live numbers.
3. `PORTING_PLAN.md` §2: the pipeline's steps and data flow, as reconstructed from the original code.

## Reference: what the code does and why

| File | What it answers |
|---|---|
| `DATA_SOURCES.md` | Where every input comes from (Zenodo, the GitHub code, reconstructed files, golden outputs). |
| `PARAMETER_TABLE.md` | The input table: rows, columns, units, what each parameter family means. |
| `PARAMETER_SELECTION.md` | How the "optimal" 41-parameter set was chosen (a 2018 optimisation, not ported). |
| `RANDOMNESS.md` | What can make two runs differ, and how the port keeps runs reproducible. |
| `PORTING_PLAN.md` | The Phase 0 audit: scope, module map, dependencies, plan. |

## Logs: running records

| File | Contents |
|---|---|
| `DECISIONS_AND_ISSUES.md` | **The master log.** Every finding and decision, by ID (A = data, B = PAD/ranking, C = groups/UMAP, D = environment, G, L, P, R, U, X; Q = open questions for the lab). Cite IDs elsewhere instead of repeating content. |
| `NUMERICAL_DIFFERENCES.md` | Every place the port's numbers differ from the original's, by how much, and the approved test tolerance (N1–N3). |
| `PERFORMANCE.md` | Speed and memory: measurements before/after each optimisation. |
| `LAB_QUESTIONS.md` | Questions for the Green/Oegema lab, written to send (the Q table in the master log tracks their status). |

## Elsewhere in the repo

| File | Contents |
|---|---|
| `../data/README.md` | The Zenodo files: record, licence, checksums. |
| `../tests/golden/README.md` | The frozen outputs of the original code, and how to regenerate them. |
| `../legacy_baseline/README.md` | How the original Python 2.7 code was made runnable unmodified. |
| `reports/` | Machine-readable Phase 1 comparisons (port/legacy vs. published files). |

## Conventions

- New finding or decision → a row in `DECISIONS_AND_ISSUES.md` (new ID). If it needs more than a few lines, write a
  topic note here and link it from the row (as B14 → `PARAMETER_SELECTION.md`).
- Mark anything inferred rather than read from code or data as **(inferred)**.
- Date entries (YYYY-MM-DD), and record who approved scientific changes.

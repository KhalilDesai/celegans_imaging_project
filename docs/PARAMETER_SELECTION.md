# How the "optimal" parameter set was chosen (and what the port does with it)

Written 2026-10-26 from a reading of the original code. Marked **(inferred)** where the code does not say so
explicitly; to be checked against the paper's STAR Methods (Green et al. 2024, PMC11166207).

## Short answer

The optimal set (17 GLS + 24 MS parameters) was chosen **once, before the published pipeline ran**, in a
separate optimisation the authors ran by hand in March 2018. The pipeline never re-runs that selection: it
reads a hardcoded result. The port reproduces everything downstream of that fixed list, not the selection itself.

## What the pipeline does (ported, `src/embd/params_use_assign.py`)

`_group_rank_opt()` (docstring dated 3/7/18) holds two hardcoded lists: `ordered_names`, about 90 candidate
parameters, and `ordered_vals`, a count from 0 to 18 for each. Parameters with a count **> 4** are removed
(`th = 4`); the remaining 17 GLS + 24 MS form the optimal set (DECISIONS B3). The docstring says the run used a
70/30 train/test split of the manual groups, at least 3 members per group in training, and that the "origin was
reset prior to this run", and it lists the training groups.

## Where the selection code lives (NOT ported, still Python 2.7)

In this repository, in the original file `development/groupAssesment.py` (unchanged since the upload of the
upstream code; there is no separate repository for it). The port's `src/embd/groupAssesment.py` lists these
functions under "NOT YET PORTED" because they are not on the path to any published output:

| Function (original file, line) | Role |
|---|---|
| `optimize_groups_by_rank` (L999) | One optimisation run. `gain` = minus the average within-group PAD rank (lower rank = group members find each other first), averaged separately over cell-fate (`cf`) and morphogenesis (`m`) groups so each half weighs equally. Optionally leaves out one group (`group_exclude`). |
| `optimizeByGain` (L823) | Greedy backward elimination: each round, try removing every remaining parameter, permanently drop the one whose removal improves the gain most; stop when no removal helps (max 150 rounds). Returns the dropped parameters. |
| `get_in_group_rank` (L1089), `calc_pad_normdist`, `get_top_hits`, `remove_param_from_use` (L791), `get_rnai_groups_obj` (L884) | Helpers for the above. |
| `evaluate_ranking` (L1109) | Scores a parameter set on the held-out test groups. |
| `optimize_groups_by_min_max` (L902), `optimizeByInteractome` (L955), `getRNAiGroupsPhenoOld` | Earlier/alternative attempts ("initial groups for first pass optimization"). |

The driver is commented out at the end of the module (around L1955):

```python
# pool = mp.Pool(processes=6)
# results_all = pool.map(optimize_groups_by_rank,
#                        zip(repeat(rlist_train), ..., [None]+range(len(rlist_train))))
```

i.e. one run with no group left out plus one run per training group left out. With 17 training groups that is
**18 runs**, which matches the maximum count of 18. So `ordered_vals` is most likely **the number of runs (out
of 18) in which a parameter was dropped**, and "count > 4" means "dropped in more than 4 of 18 runs"
**(inferred)**. How the run outputs were turned into the hardcoded lists is not in the code (presumably by hand
from the logs).

## Consequences

- **Circularity (by design).** The same hand-curated groups choose the parameters and later produce each gene's
  auto-description (group centroids). The elt-1 duplicate in group `dorsRupt` (DECISIONS A17) slightly biases both.
- **Not reproducible from the repo as is.** The optimisation is unported, its driver is commented out, the
  `[None]+range(...)` line only works in Python 2, and it ran on the 2018 state of the per-gene table, whereas
  the public data is the 2023 sheet. A re-run today would likely give a somewhat different list.
- **Porting it is a separate task**, under the open decision on non-published experiment code (with tSNE and
  the classifiers). It would allow re-deriving or testing the optimal set, but its result would not be expected to
  match the 2018 list exactly.

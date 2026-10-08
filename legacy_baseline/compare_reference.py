"""Compare legacy-baseline PAD rankings with the published all_best_matches_05192022.xlsx.

Python 3 (pandas + openpyxl). Usage:
    python compare_reference.py <data_dir> <legacy_work_Z_dir> <out_json> [<out_md>]

For each parameter set (optimal, GLS, MS) it reports:
  * exact agreement of the ordered top-10 hit list per seed,
  * set overlap of the top-10 lists and position-wise agreement,
  * |delta PAD| for every published (seed, hit) pair, looked up in the legacy all-pairs CSV,
  * |delta PAD| on the rank-aligned top-10 values.
"""
import csv
import json
import sys

import numpy as np
import pandas as pd

SETS = [('opt', 3, 5, 'phenoDist_rank_04142022.csv'),
        ('GLS', 6, 8, 'phenoDist_rank_GLS_ONLY05102022.csv'),
        ('MS', 9, 11, 'phenoDist_rank_MS_ONLY05102022.csv')]


def load_reference(path):
    df = pd.read_excel(path, header=None, dtype=object)
    df = df[pd.to_numeric(df[0], errors='coerce').notna()].copy()
    for c in range(12):
        if c not in (2, 4, 7, 10):
            df[c] = pd.to_numeric(df[c])
    return df.astype({0: int, 1: int, 3: int, 6: int, 9: int})


def load_legacy_best_matches(path):
    """Parse the legacy 'all_best_matches_<date>.csv' (repeated header lines, blank separators)."""
    rows = []
    with open(path, newline='', encoding='latin-1') as f:
        for r in csv.reader(f):
            if len(r) >= 12 and r[0].strip().isdigit():
                rows.append(r[:12])
    df = pd.DataFrame(rows)
    for c in (0, 1, 3, 6, 9):
        df[c] = df[c].astype(int)
    for c in (5, 8, 11):
        df[c] = df[c].astype(float)
    return df


def load_all_pairs(path):
    """phenoDist_rank CSV: 2 'params used' lines + header, then r1,r2,dist,dist2WT,dims,CSI,PAD."""
    df = pd.read_csv(path, skiprows=2, float_precision='round_trip')
    df['s'] = df['RNA1'].str[-4:].astype(int)
    df['h'] = df['RNA2'].str[-4:].astype(int)
    return df.set_index(['s', 'h'])['PAD']


def lists(df, idcol, padcol):
    out = {}
    for seed, g in df.groupby(1):
        g = g.sort_values(0)
        out[seed] = (g[idcol].tolist(), g[padcol].tolist())
    return out


def compare(ref, leg, pairs):
    r_l, l_l = ref, leg
    seeds = sorted(r_l)
    exact = sum(r_l[s][0] == l_l[s][0] for s in seeds)
    overlap = [len(set(r_l[s][0]) & set(l_l[s][0])) / 10.0 for s in seeds]
    pos = np.array([[a == b for a, b in zip(r_l[s][0], l_l[s][0])] for s in seeds])
    first_mismatch = [next((k + 1 for k, (a, b) in enumerate(zip(r_l[s][0], l_l[s][0])) if a != b), None)
                      for s in seeds]
    d_aligned = np.array([abs(a - b) for s in seeds for a, b in zip(r_l[s][1], l_l[s][1])])
    d_pairs = []
    missing = 0
    for s in seeds:
        for h, p in zip(*r_l[s]):
            if (s, h) in pairs.index:
                d_pairs.append(abs(float(pairs.loc[(s, h)]) - p))
            else:
                missing += 1
    d_pairs = np.array(d_pairs)
    mism = [s for s, fm in zip(seeds, first_mismatch) if fm is not None]
    return {
        'n_seeds': len(seeds),
        'exact_ordered_top10': exact,
        'exact_ordered_top10_frac': exact / float(len(seeds)),
        'mean_top10_set_overlap': float(np.mean(overlap)),
        'positionwise_agreement_by_rank': [float(x) for x in pos.mean(axis=0)],
        'seeds_with_mismatch': mism[:50],
        'first_mismatch_rank_hist': pd.Series([fm for fm in first_mismatch if fm]).value_counts().sort_index().to_dict(),
        'abs_dPAD_published_pairs': {'n': int(d_pairs.size), 'missing_pairs': missing,
                                     'max': float(d_pairs.max()), 'median': float(np.median(d_pairs)),
                                     'p99': float(np.percentile(d_pairs, 99))},
        'abs_dPAD_rank_aligned': {'max': float(d_aligned.max()), 'median': float(np.median(d_aligned))},
    }


def main(argv):
    data_dir, zdir, out_json = argv[1:4]
    out_md = argv[4] if len(argv) > 4 else None
    ref = load_reference(data_dir + '/all_best_matches_05192022.xlsx')
    leg = load_legacy_best_matches(zdir + '/all_best_matches_legacy.csv')
    res = {}
    for name, idc, padc, fn in SETS:
        pairs = load_all_pairs(zdir + '/' + fn)
        res[name] = compare(lists(ref, idc, padc), lists(leg, idc, padc), pairs)
    with open(out_json, 'w') as f:
        json.dump(res, f, indent=2, sort_keys=True, default=str)
    lines = ['| Parameter set | Exact ordered top-10 | Mean top-10 overlap | max abs dPAD | median abs dPAD | p99 abs dPAD |',
             '|---|---|---|---|---|---|']
    for name in res:
        r = res[name]
        a = r['abs_dPAD_published_pairs']
        lines.append('| %s | %d/%d (%.1f%%) | %.4f | %.3g | %.3g | %.3g |' % (
            name, r['exact_ordered_top10'], r['n_seeds'], 100 * r['exact_ordered_top10_frac'],
            r['mean_top10_set_overlap'], a['max'], a['median'], a['p99']))
    print('\n'.join(lines))
    if out_md:
        with open(out_md, 'w') as f:
            f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main(sys.argv)

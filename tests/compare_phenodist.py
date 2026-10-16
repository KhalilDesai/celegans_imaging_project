"""Compare a ported phenoDist_rank CSV with the golden legacy CSV (diagnostics, used by tests and docs).

python tests/compare_phenodist.py <golden .csv.xz> <port .csv>
"""
import json
import lzma
import sys

import numpy as np


def load(path):
    opener = lzma.open if path.endswith('.xz') else open
    with opener(path, 'rt') as f:
        lines = f.read().splitlines()
    return lines


def compare(golden_path, port_path):
    g, p = load(golden_path), load(port_path)
    res = {'n_lines': [len(g), len(p)], 'byte_identical': g == p}
    res['header_identical'] = g[:3] == p[:3]
    G = [r.split(',') for r in g[3:]]
    P = [r.split(',') for r in p[3:]]
    res['row_order_identical'] = [(a[0], a[1]) for a in G] == [(b[0], b[1]) for b in P]
    pk = {(b[0], b[1]): b for b in P}
    cols = ['dist', 'dist2WT', 'dims', 'CSI', 'PAD']
    diff_lines, maxabs, maxrel, n_text = 0, {c: 0.0 for c in cols}, {c: 0.0 for c in cols}, {c: 0 for c in cols}
    max12 = {c: 0.0 for c in cols}
    for a in G:
        b = pk[(a[0], a[1])]
        if a != b:
            diff_lines += 1
        for k, c in enumerate(cols, start=2):
            if a[k] != b[k]:
                n_text[c] += 1
                x, y = float(a[k]), float(b[k])
                if np.isfinite(x) and np.isfinite(y):
                    maxabs[c] = max(maxabs[c], abs(x - y))
                    # size of the difference in units of the 12th significant digit of the legacy value
                    unit = 10.0 ** (np.floor(np.log10(abs(x))) - 11) if x != 0 else 1e-300
                    max12[c] = max(max12[c], abs(x - y) / unit)
                    if x != 0:
                        maxrel[c] = max(maxrel[c], abs(x - y) / abs(x))
    res.update({'rows_differing': diff_lines, 'fields_differing_text': n_text,
                'max_abs_diff': maxabs, 'max_rel_diff': maxrel,
                'max_diff_in_12th_digit_units': max12})
    # ranking check: PADFilter semantics (argsort of file-order PADs, reversed, skip first), top 51
    def ranks(rows):
        by = {}
        for r in rows:
            by.setdefault(r[0], []).append(r)
        out = {}
        for s, rr in by.items():
            pads = [float(r[6]) for r in rr]
            idx = np.argsort(pads)[::-1][1:52]
            out[s] = [rr[i][1] for i in idx]
        return out
    rg, rp = ranks(G), ranks(P)
    res['seeds_top10_identical'] = sum(rg[s][:10] == rp[s][:10] for s in rg)
    res['seeds_top51_identical'] = sum(rg[s] == rp[s] for s in rg)
    res['n_seeds'] = len(rg)
    return res


if __name__ == '__main__':
    print(json.dumps(compare(sys.argv[1], sys.argv[2]), indent=1, sort_keys=True))

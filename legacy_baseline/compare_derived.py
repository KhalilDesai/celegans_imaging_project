"""Compare legacy derived columns (top_hits_PAD, reciprocal_matches, auto_description) with the
published EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv.

Python 3. Usage: python compare_derived.py <data_dir> <legacy_Z_dir> <out_json>

Gene names are mapped to EMBD IDs via the 2023 index. The published CSV uses 2022 names and has been
through Excel (DECISIONS A12), so names that don't resolve get an alias, learned from positions where
the two lists otherwise agree. The alias table is reported so it can be checked by hand.
"""
import json
import sys

import pandas as pd

WT_NOTE = 'Embryos do not exhibit an arrest phenotype'


def split(v):
    if not isinstance(v, str) or v.strip() == '':
        return []
    if v.strip() == WT_NOTE:
        return [WT_NOTE]
    return [t.strip() for t in v.split(',')]


def main(data_dir, zdir, out_json):
    pub = pd.read_csv(data_dir + '/EMBD0001_0504_updated_topmatch_recip_autodescriptions_052422.csv',
                      dtype=str, encoding='latin-1', keep_default_na=False)
    leg = pd.read_csv(zdir + '/genes_derived_legacy.csv', dtype=str, keep_default_na=False)
    idx = pd.read_excel(data_dir + '/Index_EMBD_ID_Gene_names.xlsx', dtype=str)
    name2id = {}
    for _, r in idx.iterrows():
        for n in (r['Gene Name'], r['C. elegans ID']):
            if isinstance(n, str) and n not in name2id:
                name2id[n] = int(r['EMBD ID'])
    # EMBD 87 is stored as an Excel date in the index (DECISIONS A16); the CSV has 'apr-1' or '1-Apr'.
    for n in ('apr-1', '1-Apr'):
        name2id[n] = 87
    m = pub.merge(leg, on='embd_id', suffixes=('_pub', '_leg'))
    m = m[m.embd_id.astype(int).between(1, 503)]
    # learn aliases for unresolved published names
    alias = {}
    for col in ('top_hits_PAD', 'reciprocal_matches'):
        for _, r in m.iterrows():
            a, b = split(r[col + '_pub']), split(r[col + '_leg'])
            if len(a) != len(b):
                continue
            unknown = [k for k, x in enumerate(a) if x not in name2id and x != WT_NOTE]
            # learn only if every other position agrees exactly
            if len(unknown) == 1 and all(x == y for k, (x, y) in enumerate(zip(a, b)) if k != unknown[0]):
                y = b[unknown[0]]
                if y in name2id:
                    alias.setdefault(a[unknown[0]], set()).add(y)
    alias_ids = {k: sorted({name2id[y] for y in v}) for k, v in alias.items()}

    def ids(tokens):
        out = []
        for t in tokens:
            if t == WT_NOTE:
                out.append(t)
            elif t in name2id:
                out.append(name2id[t])
            elif t in alias_ids and len(alias_ids[t]) == 1:  # unambiguous alias only
                out.append(alias_ids[t][0])
            else:
                out.append('?' + t)
        return out
    res = {'alias_table': {k: sorted(v) for k, v in alias.items()}}
    for col in ('top_hits_PAD', 'reciprocal_matches'):
        diff = []
        for _, r in m.iterrows():
            a, b = ids(split(r[col + '_pub'])), ids(split(r[col + '_leg']))
            if a != b:
                diff.append({'embd_id': int(r.embd_id), 'published': a, 'legacy': b})
        res[col] = {'n': len(m), 'identical': len(m) - len(diff), 'differing': diff}
    ad = m[m.auto_description_pub.str.strip() != m.auto_description_leg.str.strip()]
    res['auto_description'] = {'n': len(m), 'identical': len(m) - len(ad),
                               'differing_embd_ids': ad.embd_id.astype(int).tolist()}
    with open(out_json, 'w') as f:
        json.dump(res, f, indent=1, default=str)
    for col in ('top_hits_PAD', 'reciprocal_matches', 'auto_description'):
        print('%-20s identical %d/%d' % (col, res[col]['identical'], res[col]['n']))
    print('aliases:', res['alias_table'])
    return res


if __name__ == '__main__':
    main(*sys.argv[1:4])

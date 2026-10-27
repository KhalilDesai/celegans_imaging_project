"""Numbers behind docs/CONTROLS_AND_NORMALISATION.md (optimal parameter set).

usage: EMBD_OUT_DIR=<run dir>/ EMBD_SQLITE=<run dir>/embd.sqlite EMBD_PARAMSET=opt \
       PYTHONPATH=src .venv/bin/python tools/assess_controls.py
<run dir> must hold a finished run (canonical/genes_table.csv, genes_derived_port.csv), e.g. from
`pipeline reproduce` or notebooks/walkthrough.ipynb (build/walkthrough). Takes a few seconds, < 0.5 GB.
"""
import contextlib
import io
import os

import numpy as np
import pandas as pd

with contextlib.redirect_stdout(io.StringIO()):  # the original modules print the parameter lists on import
    from embd import _rnai_constants as K
    from embd.RNAiClass import PARAM_NAMES_USE_GLS as G, PARAM_NAMES_USE_MS as M, RNAiClass

OUT = os.environ['EMBD_OUT_DIR']
keys = [p + '_GLS' for p in G] + [p + '_MS' for p in M]

# 1. control embryo-to-embryo spread vs the normalisation scale (spread across RNAi conditions)
ratio = np.array([K.controlStd[k] / K.rnaiStd[k] for k in keys])
print('controlStd/rnaiStd over %d params: median %.2f, min %.2f, max %.2f, n>0.5: %d, n>1: %d'
      % (len(keys), np.median(ratio), ratio.min(), ratio.max(), (ratio > 0.5).sum(), (ratio > 1).sum()))
print('highest ratios:', [(k, round(r, 2)) for r, k in sorted(zip(ratio, keys))[-5:]])

# 2. the control row of the published table vs the hardcoded origin
t = pd.read_csv(OUT + 'canonical/genes_table.csv', low_memory=False)
c = t[t.embd_id == 0].iloc[0]
diff = []
for k in keys:
    s = k.rsplit('_', 1)[1]
    col = k if k in t.columns else k + '_' + s  # e.g. CoM0Rtd_GLS -> column CoM0Rtd_GLS_GLS
    diff.append((float(c[col]) - K.controlOrigin[k]) / K.rnaiStd[k])
print('control row vs hardcoded origin, in RNAi-spread units: max |diff| %.3g' % np.nanmax(np.abs(diff)))

# 3. where genes sit relative to the origin, by the pipeline's own WT-like / phenotype call
with contextlib.redirect_stdout(io.StringIO()):
    R = [RNAiClass(i, verb=False) for i in range(1, 504)]
    for r in R:
        r.setNDPosition()
X = np.array([r.posND for r in R])
X = X[:, ~np.all(np.isnan(X), 0)]
d2z = np.array([r.getDist2Zero() for r in R])
der = pd.read_csv(OUT + 'genes_derived_port.csv').set_index('embd_id')
call = np.array([der.loc[i, 'WT_like_or_phenotype'] for i in range(1, 504)])
for lab in ('WT-like', 'Dev. Delay', 'Phenotype'):
    m = call == lab
    print('%-10s n=%3d  distance to origin: median %.2f' % (lab, m.sum(), np.median(d2z[m])))
wt = X[call == 'WT-like']
mu = np.nanmean(wt, 0)
print('WT-like genes, per-parameter mean position: median |mean| %.2f, max %.2f, n |mean|>0.5: %d'
      % (np.median(np.abs(mu)), np.abs(mu).max(), (np.abs(mu) > 0.5).sum()))
print('WT-like genes, per-parameter std: median %.2f' % np.median(np.nanstd(wt, 0)))

# 4. redundancy between parameters
C = pd.DataFrame(X).corr().values
iu = np.triu_indices(C.shape[0], 1)
print('pairwise |r|: median %.2f, pairs |r|>0.8: %d of %d'
      % (np.median(np.abs(C[iu])), (np.abs(C[iu]) > 0.8).sum(), len(iu[0])))
ev = np.linalg.svd(np.nan_to_num(X - np.nanmean(X, 0)), compute_uv=False) ** 2
ev /= ev.sum()
print('PCA variance in first 1/3/5/10 components: %.2f %.2f %.2f %.2f; participation ratio %.1f'
      % (ev[:1].sum(), ev[:3].sum(), ev[:5].sum(), ev[:10].sum(), 1 / np.sum(ev ** 2)))

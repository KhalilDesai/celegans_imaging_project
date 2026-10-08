"""Run the ORIGINAL development/makeUMAPs_V2.py (unmodified) and record the UMAP embedding.

The script runs at module level and only saves an SVG. To capture coordinates, ``umap.UMAP`` is
wrapped so that ``fit_transform`` saves its output. The call itself (default arguments, global
``np.random.seed(71076)`` set by the script) is unchanged. Inputs the script reads are expected in
<work>/Z:/: ``gene_key.csv`` and ``PAD_Matrix_042122_optParams.npy`` (copied from the legacy
matrices step by the caller, LEGACY-DRIVER).

Usage: python run_umap_legacy.py <tag>      # tag names the output files, e.g. opt
"""
import json
import os
import runpy
import sys

import numpy as np
import umap

REPO = os.environ.get('EMBD_REPO', '/repo')
WORK = os.environ.get('EMBD_WORK', '/work')


def main(tag):
    os.chdir(WORK)
    captured = []
    Orig = umap.UMAP

    class RecordingUMAP(Orig):
        def fit_transform(self, X, *a, **k):
            out = super(RecordingUMAP, self).fit_transform(X, *a, **k)
            captured.append((self.get_params(), np.asarray(out)))
            return out
    umap.UMAP = RecordingUMAP
    runpy.run_path(os.path.join(REPO, 'development', 'makeUMAPs_V2.py'), run_name='__main__')
    params, pts = captured[-1]
    np.save(os.path.join('Z:', 'umap_points_%s.npy' % tag), pts)
    with open(os.path.join('Z:', 'umap_params_%s.json' % tag), 'w') as f:
        json.dump({k: (v if isinstance(v, (int, float, str, bool, type(None))) else repr(v))
                   for k, v in params.items()}, f, indent=1, sort_keys=True)
    with open('/legacy-umap-versions.txt') as f:
        versions = f.read()
    with open(os.path.join('Z:', 'umap_versions_%s.txt' % tag), 'w') as f:
        f.write(versions)
    print('saved', pts.shape)


if __name__ == '__main__':
    main(sys.argv[1])

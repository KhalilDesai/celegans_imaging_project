"""Freeze the legacy-baseline outputs as golden reference files (DECISIONS X4, option (a)).

Usage: python tests/golden/build_golden.py <legacy_Z_dir> <umap_work_dir> <out_dir>

Large files are stored xz-compressed (exact bytes preserved). Writes:
  MANIFEST.sha256    sha256 of every stored file (verified by the test suite)
  CONTENT.sha256     sha256 of the *uncompressed* content (what a re-run must reproduce)
  provenance.json    how and from what the goldens were produced
"""
import glob
import hashlib
import json
import lzma
import os
import shutil
import subprocess
import sys

SETS = ('opt', 'gls', 'ms', 'all')
PHENODIST = {'opt': 'phenoDist_rank_04142022.csv', 'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
             'ms': 'phenoDist_rank_MS_ONLY05102022.csv', 'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def main(zdir, umapdir, out):
    os.makedirs(out, exist_ok=True)
    plan = []  # (source, dest_name, compress)
    for s in SETS:
        plan += [(os.path.join(zdir, PHENODIST[s]), 'phenodist_%s.csv' % s, True),
                 (os.path.join(zdir, 'PAD_Matrix_%s.npy' % s), 'pad_matrix_%s.npy' % s, True),
                 (os.path.join(zdir, 'Dist_Matrix_%s.npy' % s), 'dist_matrix_%s.npy' % s, True),
                 (os.path.join(zdir, 'posND_%s.npy' % s), 'posND_%s.npy' % s, True),
                 (os.path.join(zdir, 'dist2zero_%s.npy' % s), 'dist2zero_%s.npy' % s, False),
                 (os.path.join(zdir, 'params_used_%s.json' % s), 'params_used_%s.json' % s, False),
                 (os.path.join(umapdir, s, 'Z:', 'umap_points_%s.npy' % s), 'umap_points_%s.npy' % s, False),
                 (os.path.join(umapdir, s, 'Z:', 'umap_params_%s.json' % s), 'umap_params_%s.json' % s, False)]
    plan += [(os.path.join(zdir, 'all_best_matches_legacy.csv'), 'all_best_matches.csv', False),
             (os.path.join(zdir, 'hits_matrix_v2_legacy.csv'), 'hits_matrix_v2.csv', False),
             (os.path.join(zdir, 'pad_rank_legacy.csv'), 'pad_rank.csv', False),
             (os.path.join(zdir, 'genes_derived_legacy.csv'), 'genes_derived.csv', False),
             (os.path.join(zdir, 'Matrix_labels_opt.npy'), 'matrix_labels.npy', False),
             (os.path.join(zdir, 'gene_key.csv'), 'gene_key.csv', False),
             (os.path.join(umapdir, 'opt', 'Z:', 'umap_versions_opt.txt'), 'umap_versions.txt', False)]
    content = {}
    for src, name, compress in plan:
        if not os.path.exists(src):
            raise SystemExit('missing legacy output: %s' % src)
        content[name] = sha(src)
        if compress:
            with open(src, 'rb') as fi, lzma.open(os.path.join(out, name + '.xz'), 'wb', preset=9) as fo:
                shutil.copyfileobj(fi, fo)
        else:
            shutil.copyfile(src, os.path.join(out, name))
    with open(os.path.join(out, 'CONTENT.sha256'), 'w') as f:
        for k in sorted(content):
            f.write('%s  %s\n' % (content[k], k))

    def docker(cmd):
        try:
            return subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
        except Exception:
            return None
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    prov = {
        'description': 'Outputs of the ORIGINAL Python 2.7 downstream code (unmodified) run by '
                       'legacy_baseline/run_legacy.py on data/All_genes_mySQL_table.xlsx. '
                       'See docs/DECISIONS_AND_ISSUES.md sections L and X.',
        'legacy_image': docker(['docker', 'image', 'inspect', 'embd-legacy:py27', '--format', '{{.Id}}']),
        'legacy_umap_image': docker(['docker', 'image', 'inspect', 'embd-legacy:umap-py36', '--format', '{{.Id}}']),
        'legacy_python_packages': docker(['docker', 'run', '--rm', '--platform', 'linux/amd64', '--entrypoint', 'cat',
                                          'embd-legacy:py27', '/legacy-versions.txt']),
        'git_commit_of_harness': docker(['git', '-C', repo, 'rev-parse', 'HEAD']),
        'input_sheet_sha256': sha(os.path.join(repo, 'data', 'All_genes_mySQL_table.xlsx')),
        'gene_key_note': 'gene_key.csv reconstructed from the 2023 Zenodo index (DECISIONS A9, A12, A16)',
        'timings_legacy': open(os.path.join(zdir, 'timings_legacy.jsonl')).read().splitlines(),
    }
    with open(os.path.join(out, 'provenance.json'), 'w') as f:
        json.dump(prov, f, indent=1, sort_keys=True)
    stored = sorted(p for p in glob.glob(os.path.join(out, '*')) if os.path.basename(p) not in
                    ('MANIFEST.sha256', 'build_golden.py', '__pycache__', 'README.md'))
    with open(os.path.join(out, 'MANIFEST.sha256'), 'w') as f:
        for p in stored:
            f.write('%s  %s\n' % (sha(p), os.path.basename(p)))
    print('stored %d files' % len(stored))


if __name__ == '__main__':
    main(*sys.argv[1:4])

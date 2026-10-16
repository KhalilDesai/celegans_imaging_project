"""Regression tests: Python 3 port vs. golden legacy outputs (tests/golden).

Tolerances and their justifications: docs/NUMERICAL_DIFFERENCES.md.

* Fast tier (default): builds the emulated DB, checks every gene's ND vector bit-for-bit and PAD
  values / top-10 rankings for a sample of seeds. About a minute.
* Slow tier (``pytest -m slow``): compares every output of a full pipeline run. Set
  ``EMBD_PORT_RUN_DIR`` to an existing ``python -m embd.pipeline all`` output directory to reuse it;
  otherwise the full pipeline is run (several minutes).
"""
import csv
import io
import json
import lzma
import os
import subprocess
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(REPO, 'tests', 'golden')
SRC = os.path.join(REPO, 'src')
sys.path.insert(0, SRC)
sys.path.insert(0, os.path.join(REPO, 'tests'))

SETS = ('opt', 'gls', 'ms', 'all')
PHENODIST = {'opt': 'phenoDist_rank_04142022.csv', 'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
             'ms': 'phenoDist_rank_MS_ONLY05102022.csv', 'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}
PAD_ATOL_TEXT = 2e-12   # N1: two units in the 12th printed digit
MATRIX_RTOL = 1e-13     # N1: full-precision PAD / distance matrices
SAMPLE_SEEDS = [1, 2, 3, 13, 45, 76, 87, 117, 181, 185, 255, 310, 363, 408, 434, 478, 503]


def golden_npy(name):
    path = os.path.join(GOLDEN, name)
    if os.path.exists(path + '.xz'):
        return np.load(io.BytesIO(lzma.open(path + '.xz').read()), allow_pickle=True)
    return np.load(path, allow_pickle=True)


def golden_lines(name):
    path = os.path.join(GOLDEN, name)
    if os.path.exists(path + '.xz'):
        return lzma.open(path + '.xz', 'rt').read().splitlines()
    with open(path) as f:
        return f.read().splitlines()


def _env(out_dir, paramset=None):
    env = dict(os.environ)
    env.update({'EMBD_OUT_DIR': str(out_dir), 'EMBD_SQLITE': os.path.join(str(out_dir), 'embd.sqlite'),
                'PYTHONPATH': SRC + os.pathsep + env.get('PYTHONPATH', '')})
    if paramset:
        env['EMBD_PARAMSET'] = paramset
    return env


def _canonical(tmp_root):
    """Canonical genes CSV + gene key, built from data/ (exactly as the legacy run used them)."""
    from embd_data import sheet
    d = os.path.join(tmp_root, 'canonical')
    os.makedirs(d, exist_ok=True)
    sheet.convert_genes_table(os.path.join(REPO, 'data', 'All_genes_mySQL_table.xlsx'), os.path.join(d, 'genes_table.csv'))
    sheet.write_gene_key(os.path.join(REPO, 'data', 'Index_EMBD_ID_Gene_names.xlsx'), os.path.join(d, 'gene_key.csv'))
    with open(os.path.join(d, 'gene_key.csv'), 'rb') as a, open(os.path.join(GOLDEN, 'gene_key.csv'), 'rb') as b:
        assert a.read() == b.read(), 'gene_key reconstruction differs from the one used for the goldens'
    return d


def _run(args, out_dir, paramset=None):
    subprocess.check_call([sys.executable, '-m', 'embd.pipeline'] + args, env=_env(out_dir, paramset),
                          stdout=subprocess.DEVNULL)


@pytest.fixture(scope='session')
def fast_run(tmp_path_factory):
    root = tmp_path_factory.mktemp('fast')
    canon = _canonical(str(root))
    out = root / 'out'
    out.mkdir()
    _run(['init-db', '--genes-csv', os.path.join(canon, 'genes_table.csv'), '--gene-key', os.path.join(canon, 'gene_key.csv')], out)
    for s in SETS:
        _run(['vectors', s], out, paramset=s)
    return out


@pytest.mark.parametrize('ps', SETS)
def test_nd_vectors_bit_identical(fast_run, ps):
    g = golden_npy('posND_%s.npy' % ps)
    p = np.load(os.path.join(str(fast_run), 'posND_%s.npy' % ps))
    assert np.array_equal(np.isnan(g), np.isnan(p))
    assert np.array_equal(g[~np.isnan(g)].view(np.int64), p[~np.isnan(p)].view(np.int64))
    assert json.load(open(os.path.join(GOLDEN, 'params_used_%s.json' % ps))) == \
        json.load(open(os.path.join(str(fast_run), 'params_used_%s.json' % ps)))


@pytest.mark.parametrize('ps', SETS)
def test_dist2zero_within_tolerance(fast_run, ps):
    g = golden_npy('dist2zero_%s.npy' % ps)
    p = np.load(os.path.join(str(fast_run), 'dist2zero_%s.npy' % ps))
    np.testing.assert_allclose(p, g, rtol=MATRIX_RTOL, atol=0)


@pytest.mark.parametrize('ps', SETS)
def test_sample_pads_and_rankings(fast_run, ps):
    """PAD for every pair involving the sample seeds, recomputed by the port's RNAiClass, vs the golden matrix;
    and PADFilter top-10 from the golden CSV text reproduced from the port's values."""
    code = r'''
import json, sys, numpy as np
from embd.RNAiClass import RNAiClass
seeds = json.loads(sys.argv[1])
objs = {}
for i in range(1, 504):
    r = RNAiClass(i, verb=False); r.setNDPosition(); objs[i] = r
P = np.array([[objs[s].getPAD(objs[j]) for j in range(1, 504)] for s in seeds])
np.save(sys.argv[2], P)
'''
    outp = os.path.join(str(fast_run), 'sample_pad_%s.npy' % ps)
    subprocess.check_call([sys.executable, '-c', code, json.dumps(SAMPLE_SEEDS), outp], env=_env(fast_run, ps),
                          stdout=subprocess.DEVNULL)
    P = np.load(outp)
    G = golden_npy('pad_matrix_%s.npy' % ps)[[s - 1 for s in SAMPLE_SEEDS]]
    np.testing.assert_allclose(P, G, rtol=MATRIX_RTOL, atol=1e-15)
    for row_p, row_g in zip(P, G):
        assert list(np.argsort(row_p)[::-1][:11]) == list(np.argsort(row_g)[::-1][:11])


# ----------------------------------------------------------------------------------------- slow
@pytest.fixture(scope='session')
def full_run(tmp_path_factory):
    existing = os.environ.get('EMBD_PORT_RUN_DIR')
    if existing:
        return existing
    root = tmp_path_factory.mktemp('full')
    canon = _canonical(str(root))
    out = root / 'out'
    out.mkdir()
    _run(['init-db', '--genes-csv', os.path.join(canon, 'genes_table.csv'), '--gene-key', os.path.join(canon, 'gene_key.csv')], out)
    _run(['all'], out)
    return str(out)


@pytest.mark.slow
@pytest.mark.parametrize('ps', SETS)
def test_phenodist_csv(full_run, ps):
    from compare_phenodist import compare
    r = compare(os.path.join(GOLDEN, 'phenodist_%s.csv.xz' % ps), os.path.join(full_run, PHENODIST[ps]))
    assert r['header_identical'] and r['row_order_identical']
    assert r['n_lines'][0] == r['n_lines'][1]
    for c in ('dist', 'dist2WT', 'dims', 'CSI'):
        assert r['fields_differing_text'][c] == 0, (c, r)
    assert r['max_abs_diff']['PAD'] <= PAD_ATOL_TEXT, r
    assert r['seeds_top51_identical'] == r['n_seeds'] == 503


@pytest.mark.slow
@pytest.mark.parametrize('ps', SETS)
@pytest.mark.parametrize('kind', ['pad_matrix', 'dist_matrix'])
def test_matrices(full_run, ps, kind):
    g = golden_npy('%s_%s.npy' % (kind, ps))
    p = np.load(os.path.join(full_run, '%s_%s.npy' % ({'pad_matrix': 'PAD_Matrix', 'dist_matrix': 'Dist_Matrix'}[kind], ps)))
    assert g.shape == p.shape
    np.testing.assert_allclose(p, g, rtol=MATRIX_RTOL, atol=1e-15)


@pytest.mark.slow
@pytest.mark.parametrize('golden_name,port_name', [
    ('all_best_matches.csv', 'all_best_matches_port.csv'),
    ('hits_matrix_v2.csv', 'hits_matrix_v2_port.csv'),
    ('pad_rank.csv', 'pad_rank_port.csv'),
    ('genes_derived.csv', 'genes_derived_port.csv'),
])
def test_derived_outputs_identical(full_run, golden_name, port_name):
    """Rankings, match lists and group assignments must match exactly."""
    g = golden_lines(golden_name)
    with open(os.path.join(full_run, port_name), encoding='utf-8') as f:
        p = f.read().splitlines()
    if golden_name in ('pad_rank.csv', 'genes_derived.csv'):
        g, p = list(csv.reader(g)), list(csv.reader(p))
    assert len(g) == len(p)
    bad = [(i, a, b) for i, (a, b) in enumerate(zip(g, p)) if a != b]
    assert not bad, bad[:5]

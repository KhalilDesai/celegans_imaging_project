"""Fast unit tests for the Python 3 port (seconds)."""
import importlib.util
import io
import lzma
import math
import os
import struct
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(REPO, 'tests', 'golden')
sys.path.insert(0, os.path.join(REPO, 'src'))


def _bits(x):
    return struct.pack('<d', x)


def test_rnai_constants_bit_identical_to_original_source():
    spec = importlib.util.spec_from_file_location('extract', os.path.join(REPO, 'tools', 'extract_rnai_constants.py'))
    ext = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ext)
    src = open(os.path.join(REPO, 'development', 'RNAiClass.py'), encoding='latin-1').read()
    original = ext.extract(src)
    from embd import _rnai_constants as K
    for name in ext.NAMES:
        port = getattr(K, name)
        assert list(port) == list(original[name]), name  # same keys, same order
        for k in port:
            a, b = port[k], original[name][k]
            assert (math.isnan(a) and math.isnan(b)) or _bits(a) == _bits(b), (name, k)


def test_param_names_identical_to_original():
    src = open(os.path.join(REPO, 'embryos_analyze', 'varLookup.py')).read()
    ns = {}
    exec(src[src.index('PARAM_NAMES = ['):src.index('# movement and alignment validation file')], ns)
    from embd.varLookup import PARAM_NAMES
    assert PARAM_NAMES == ns['PARAM_NAMES']


@pytest.mark.parametrize('ps,n_gls,n_ms', [('opt', 17, 24), ('gls', 58, 0), ('ms', 0, 50), ('all', 58, 50)])
def test_parameter_sets(monkeypatch, ps, n_gls, n_ms):
    monkeypatch.setenv('EMBD_PARAMSET', ps)
    from embd import params_use_assign
    g, m = params_use_assign.get_params_use()
    assert (len(g), len(m)) == (n_gls, n_ms)


def test_py2_str_matches_golden_csv_text():
    """Every numeric field of the legacy CSV must be reproduced by py2_str(float(field))."""
    from embd._compat import py2_str
    with lzma.open(os.path.join(GOLDEN, 'phenodist_opt.csv.xz'), 'rt') as f:
        lines = f.read().splitlines()[3:]
    n = 0
    for line in lines[::97]:  # a 2.6k-row sample; the slow test checks every byte
        fields = line.split(',')
        for s in fields[2:]:
            v = int(s) if s.lstrip('-').isdigit() else float(s)
            assert py2_str(v) == s, s
            n += 1
    assert n > 10000


def test_py2_str_cases():
    from embd._compat import py2_str
    assert py2_str(1.0) == '1.0'
    assert py2_str(0) == '0'
    assert py2_str(0.0) == '0.0'
    assert py2_str(np.float64(0.811731716522123)) == '0.811731716522'
    assert py2_str(float('inf')) == 'inf'
    assert py2_str(1e-05) == '1e-05'
    assert py2_str(123456789012345.0) == '1.23456789012e+14'


def test_getNDDistance_basic():
    from embd.embdFunc import getNDDistance
    from embd.varLookup import PARAM_NAMES
    n = len(PARAM_NAMES)
    v1 = np.full(2 * n, np.nan)
    v2 = np.full(2 * n, np.nan)
    use = np.array([['maxR', 'maxG'], ['aG']], dtype=object)
    v1[0], v1[1], v1[n + 2] = 3., 4., 1.
    v2[0], v2[1] = 0., 0.
    # overlap = 2 coords; v1 has 1 extra -> penalty sqrt(3/(3-1-0))
    assert getNDDistance(v1, v2, use) == pytest.approx(np.sqrt(3. / 2.) * 5.)
    assert np.isinf(getNDDistance(v1, np.full(2 * n, np.nan), use))


def test_out_of_scope_paths_raise():
    from embd import RNAiClass as R
    from embd import myFunc
    with pytest.raises(NotImplementedError, match='image analysis not yet ported'):
        myFunc.loadImTif('x.tif')
    with pytest.raises(NotImplementedError, match='image analysis not yet ported'):
        R.get_control_origin()


def test_padfilter_row_index_equals_original_np_where():
    """The output-neutral label index must select exactly the rows the original np.where selected."""
    from embd import PADassist
    rng = np.random.default_rng(0)
    data = [['EMBD%04d' % rng.integers(1, 6), 'EMBD%04d' % j, 1.0, 2.0, 3.0, 0.0, float(rng.random())]
            for j in range(200)]
    for name in ['EMBD0001', 'EMBD0003', 'EMBD0005', 'EMBD0009']:
        expected = np.where(np.array(data)[:, 0] == name)[0]
        assert np.array_equal(PADassist._rows_for(data, name), expected)
        orig_inds = np.where(np.array(data)[:, 0] == name)[0]
        dRNA = [data[i] for i in orig_inds]
        pads = [p[6] for p in dRNA]
        orig = [(dRNA[i][1], dRNA[i][6], dRNA[i][4]) for i in np.argsort(pads)[::-1][1:4]]
        assert PADassist.PADFilter(name, data, 3) == orig


def test_umap_seeded_is_deterministic():
    """N3: with the default random_state the port's UMAP is bit-reproducible."""
    pytest.importorskip('umap')
    from embd import makeUMAPs_V2
    rng = np.random.default_rng(1)
    X = rng.random((120, 30))
    # without svg_path, run() only computes the embedding (the 503-gene group overlay is for plotting)
    e1 = makeUMAPs_V2.run(X, query=1)
    e2 = makeUMAPs_V2.run(X, query=1)
    assert np.array_equal(e1, e2)

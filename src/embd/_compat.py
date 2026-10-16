"""Helpers that reproduce Python 2.7 / numpy 1.13 behaviour where it reaches the outputs."""
import math

NOT_PORTED = "image analysis not yet ported"


def not_ported(name):
    def _f(*args, **kwargs):
        raise NotImplementedError("%s (%s)" % (NOT_PORTED, name))
    _f.__name__ = name
    _f.__qualname__ = name
    return _f


def py2_str(x):
    """``str(x)`` / ``'{}'.format(x)`` as Python 2.7 + numpy 1.13 render numbers.

    Python 2 ``str(float)`` is ``'%.12g'`` and keeps a trailing ``'.0'`` for integral values.
    numpy 1.13 ``float64.__str__`` does the same. Ints, bools and strings pass through.
    The original CSV outputs (phenoDist_rank_*.csv, all_best_matches_*.csv) were written this
    way, so rankings were computed from 12-significant-digit values (DECISIONS B6, A11).
    """
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, int):
        return str(x)
    if isinstance(x, str):
        return x
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    if hasattr(x, 'dtype') and getattr(x.dtype, 'kind', 'f') in 'iu':
        return str(int(x))
    if math.isnan(f):
        return 'nan'
    if math.isinf(f):
        return 'inf' if f > 0 else '-inf'
    s = '%.12g' % f
    if '.' not in s and 'e' not in s and 'n' not in s:
        s += '.0'
    return s

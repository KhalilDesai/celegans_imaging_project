"""Port of embryos_analyze/myMath.py: the helpers used by the downstream code.

Curve fitting / geometry used only by per-embryo analysis (fitLine, fitSigm, fitSigmBump,
calcMoment, getCenterOfMass, AvgCurve) is out of scope and raises NotImplementedError.
"""
import numpy as np

from ._compat import not_ported

fitLine = not_ported('myMath.fitLine')
fitSigm = not_ported('myMath.fitSigm')
fitSigmBump = not_ported('myMath.fitSigmBump')
calcMoment = not_ported('myMath.calcMoment')
getCenterOfMass = not_ported('myMath.getCenterOfMass')
AvgCurve = not_ported('myMath.AvgCurve')


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def getSubIndex(a, left, right):
    if left > right: left, right = right, left
    if a[0] > a[-1]:
        return a.size - getSubIndex(a[::-1], left, right) - 1
    else:
        return np.where(a[np.where(a <= right)[0]] >= left)[0]


def getSig(x):
    x = abs(x)
    if x != np.inf and not np.isnan(x) and x != 0:
        return -int(np.floor(np.log10(x)))
    else:
        return 0

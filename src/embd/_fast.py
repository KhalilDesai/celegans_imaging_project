"""Row-vectorised pairwise distance / PAD / overlap, bit-identical to the original per-pair methods.

PERFORMANCE #4. For gene i against every gene j this returns exactly what
``RNAi[i].getDistance(RNAi[j])``, ``RNAi[i].getPAD(RNAi[j])`` and ``RNAi[i].getOverlapDims(RNAi[j])``
return, including value types (``np.inf`` / ``0.`` / int ``0``, which print differently in the CSV).

Why bit-identical: the masks, overlap counts, penalty factor, division and ``1 - ratio`` are
element-wise IEEE operations (correctly rounded, so array and scalar evaluation agree). The only
order-sensitive operation, the sum of squares in the Euclidean norm, is still done per pair with the
same BLAS call that ``np.linalg.norm`` uses for a 1-D vector (``sqrt(x.dot(x))``). Verified by
byte-comparing the phenoDist CSV and the PAD/distance matrices with the previous implementation.
"""
import numpy as np


class PairwiseRows(object):
    def __init__(self, RNAi):
        from . import RNAiClass as R
        for r in RNAi:
            if not r.NDFlag:
                r.setNDPosition()
        self.RNAi = RNAi
        self.P = np.array([r.posND for r in RNAi])
        self.nonnan_count = np.invert(np.isnan(self.P)).sum(axis=1)
        self.d20 = [r.getDist2Zero() for r in RNAi]
        # getOverlapDims counts names in both filtered paramsUse lists; those lists have unique entries
        # (dict keys), so the count equals the size of the set intersection.
        self.used = [frozenset([('GLS', pN) for pN in r.paramsUseGLS if pN in R.PARAM_NAMES_USE_GLS] +
                               [('MS', pN) for pN in r.paramsUseMS if pN in R.PARAM_NAMES_USE_MS]) for r in RNAi]

    def row(self, i):
        """(dist, pad, overlap) lists for gene i against every gene j, in j order."""
        r = self.RNAi[i]
        v1 = self.P[i]
        N = len(r.pNamesUse[0]) + len(r.pNamesUse[1])
        ov = np.invert(np.isnan(v1 * self.P))
        ovc = ov.sum(axis=1)
        i1 = self.nonnan_count[i] - ovc
        i2 = self.nonnan_count - ovc
        diff = v1 - self.P
        d20_i = self.d20[i]
        dists, pads, overlaps = [], [], []
        for j in range(len(self.RNAi)):
            if ovc[j] == 0:
                dist = np.inf
            else:
                x = diff[j][ov[j]]
                dist = np.sqrt(1. * N / (N - i1[j] - i2[j])) * np.sqrt(x.dot(x))
            if np.isinf(dist):
                pad = 0.
            else:
                with np.errstate(invalid='ignore', divide='ignore'):
                    pad = 1. - dist / (d20_i + self.d20[j])
                if not pad > 0:
                    pad = 0
            dists.append(dist)
            pads.append(pad)
            overlaps.append(len(self.used[i] & self.used[j]))
        return dists, pads, overlaps


def getPADS_all(RNAi):
    """Same list as ``[phenoDistCalc.getPADS((RNAi, i)) for i in range(len(RNAi))]``."""
    rows = PairwiseRows(RNAi)
    distAll = []
    for i in range(len(RNAi)):
        print('distance', i)
        dists, pads, overlaps = rows.row(i)
        distRNAi = [(RNAi[i].label, RNAi[j].label, dists[j], rows.d20[i], overlaps[j], pads[j])
                    for j in range(len(RNAi))]
        distRNAi.sort(key=lambda x: x[2])
        distAll.append(distRNAi[::-1])
    return distAll

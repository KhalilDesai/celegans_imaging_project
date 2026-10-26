"""Profile the slow pipeline steps in a single process (no multiprocessing, so cProfile sees everything).

Usage: EMBD_OUT_DIR=... EMBD_SQLITE=... python tools/profile_steps.py <step> [n_seeds]
steps: pads (getPADS for n seeds), csi (getCSIS for n seeds), bestmatches (n seeds), matrices (n rows)
"""
import cProfile
import io
import os
import pstats
import sys
import time

step, n = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20


def objs():
    from embd.phenoDistCalc import get_RNAi_object_list
    return get_RNAi_object_list()


def run():
    if step == 'pads':
        from embd import phenoDistCalc as P
        R = objs()
        t = time.time()
        for i in range(n):
            P.getPADS((R, i))
        return time.time() - t, 503
    if step == 'csi':
        from embd import phenoDistCalc as P
        R = objs()
        distAll = [P.getPADS((R, i)) for i in range(len(R))]
        t = time.time()
        for i in range(n):
            P.getCSIS((distAll, i))
        return time.time() - t, 503
    if step == 'bestmatches':
        from embd import PADassist
        PADassist._load = None
        t = time.time()
        PADassist.output_all_best_matches_to_csv_combined(os.environ['EMBD_OUT_DIR'] + '/prof_abm.csv', range(1, n + 1))
        return time.time() - t, 503
    if step == 'matrices':
        R = objs()
        t = time.time()
        for i in range(n):
            [R[i].getPAD(R[j]) for j in range(len(R))]
            [R[i].getDistance(R[j]) for j in range(len(R))]
        return time.time() - t, 503


pr = cProfile.Profile()
pr.enable()
dt, total = run()
pr.disable()
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(18)
print(s.getvalue()[:6000])
print('STEP %s: %.2f s for %d units -> extrapolated %.1f s for %d (single core)' % (step, dt, n, dt / n * total, total))

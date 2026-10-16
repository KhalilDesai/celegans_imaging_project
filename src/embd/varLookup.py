"""Port of embryos_analyze/varLookup.py (in-scope parts: PARAM_NAMES, constants, printLog)."""
import json
import os
import time

from . import config

nCores = 6  # number of parallel processes to run
FOLDER_IN = config.folder_in()

PARAM_NAMES = ['maxR', 'maxG', 'aG', 'aGErr', 'bG', 'bGErr', 'mG', 'mGErr', 'sG', 'sGErr', 'rG', 'rGErr', 'Gtail',
               'GtailErr', 'residsEarlyG', 'residsLateG',
               'aR', 'aRErr', 'bR', 'bRErr', 'mR', 'mRErr', 'sR', 'sRErr', 'rR', 'rRErr', 'Rtail', 'RtailErr',
               'residsEarlyR', 'residsLateR',
               'aY', 'aYErr', 'bY', 'bYErr', 'mY', 'mYErr', 'sY', 'sYErr', 'rY', 'rYErr', 'Ytail', 'YtailErr',
               'residsEarlyY', 'residsLateY',
               'tMove', 'movement', 'mRmG', 'mYmG', 'mYmR', 'tScale'] + \
              ['fracR', 'fracG', 'fracY'] + \
              ['aSigHead', 'mSigHead', 'sSigHead', 'aGaussHead', 'mGaussHead', 'sGaussHead', 'maxHead', 'tailHead',
               'scaleHead', 'devHead'] + \
              ['scaleLength', 'devLength', 'tailLength'] + \
              ['CoM{0}{1}{2}_GLS'.format(ax, c, par) for c in ['R', 'Y'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['MoI{0}{1}{2}_GLS'.format(ax, c, par) for c in ['G', 'R', 'Y'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['CoM{0}{1}{2}_MS'.format(ax, c, par) for c in ['G', 'R'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']] + \
              ['MoI{0}{1}{2}_MS'.format(ax, c, par) for c in ['G', 'R'] for ax in range(2) for par in
               ['td', 'scale', 'avgRes', 'stdTail']]

# Constants used by the per-embryo analysis (kept for completeness; only DT0_GLS_2_MS is
# imported by in-scope modules).
RADIUS_AVG = 53.
LENGTH_AVG = 100.
AG_AVG_GLS = 90
AR_AVG_GLS = 110
AY_AVG_GLS = 82
SG_AVG_GLS = 1.3
SR_AVG_GLS = 1.4
SY_AVG_GLS = 1.4
T0R_M_T0G_AVG_GLS = 3.0
T0Y_M_T0G_AVG_GLS = 3.7
T0Y_M_T0R_AVG_GLS = 0.7
GLS_MOVE_PRECISION = 4
T0_2_MOVE_AVG_GLS = 17
AR_AVG_MS = 6.1e8
AG_AVG_MS = 4.4e8
SR_AVG_MS = 2.6
SG_AVG_MS = 2.2
AVG_R_BASE_MS = 8.7e7
AVG_G_BASE_MS = 2.1e8
T0G_M_T0R_AVG_MS = 3
T0_2_MOVE_AVG_MS = 12
MS_MOVE_PRECISION = 4
MG_2_MOVE_AVG_MS = 8
MG_2_MOVE_STD_MS = 2
DT0_GLS_2_MS = T0_2_MOVE_AVG_GLS - T0_2_MOVE_AVG_MS
nLambda_MS = 6.
tauScale_MS = 0.0
nLambda_GLS = 3.
tauScale_GLS = 0.


def get_log_file():
    return FOLDER_IN + 'LOG/logFile{0}.txt'.format(time.strftime('%m_%d_%Y'))


def printLog(string):
    print(string)
    try:
        os.makedirs(FOLDER_IN + 'LOG', exist_ok=True)
        with open(get_log_file(), "a") as myfile:
            myfile.write(time.strftime("%a, %d %b %Y %H:%M:%S || ", time.localtime()) + json.dumps(string) + '\n')
    except Exception as e:
        print('Can not print to log file Error:{0}'.format(str(e)))

"""Port of embryos_analyze/myFunc.py: the image-free helpers used by the downstream code.

Image I/O and processing functions (loadImFolder, loadImTif, maxIntensProject, saveImages,
saveImagesMulti, getLargestObj, test-image generators) are out of scope and raise
NotImplementedError. The originals import cv2/tifffile/PIL at module level; this module does not.
"""
import errno
import os
import re

import numpy as np

from ._compat import not_ported

loadImFolder = not_ported('myFunc.loadImFolder')
loadImTif = not_ported('myFunc.loadImTif')
loadTestSphereIm = not_ported('myFunc.loadTestSphereIm')
loadTestCylinderIm = not_ported('myFunc.loadTestCylinderIm')
loadTestCylinderImEndOn = not_ported('myFunc.loadTestCylinderImEndOn')
maxIntensProject = not_ported('myFunc.maxIntensProject')
saveImages = not_ported('myFunc.saveImages')
saveImagesMulti = not_ported('myFunc.saveImagesMulti')
getLargestObj = not_ported('myFunc.getLargestObj')
clearFolder = not_ported('myFunc.clearFolder')


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tryint(s):
    try:
        return int(s)
    except:  # noqa: E722 (original)
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def sort_nicely(l):  # noqa: E741 (original name)
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)


def find_consec_inds(inds, n):
    """ finds consecutive indicies in a list
    INPUT:
    inds: a list of indicies
    n: the number of consecutive indicies to check for
    OUTPUT:
    returns True if consecutive indices are found and returns false if they are not. Returns pos: index position(s) where n consecutive inds are present
    """

    all_rolls = []
    a = np.roll(inds, 1)  # first roll by one
    b = np.abs(inds - a)  # subtract to find which indicies differ by 1
    c = (b == 1)  # converts to  true, false list (0,1,1,0,1)
    all_rolls.append(c)  # appends the first roll
    for i in range(n - 2):
        d = np.roll(c, 1)
        all_rolls.append(d)
        c = d
    sum_all = np.sum(all_rolls, axis=0)
    if n-1 not in sum_all:
        return False, []
    else:
        end_pos = np.where(sum_all >= n-1)[0]
        pos = end_pos - (n-1)
        return True, pos


def readCSV(fileName, delimiter=' '):
    import csv

    with open(fileName, 'r', newline='') as f:  # py2 original: csv.reader(open(fileName, 'r'))
        csvFile = csv.reader(f, delimiter=delimiter)
        result = []
        for row in csvFile:
            result.append([])
            for val in row:
                try:
                    if val != 'None' and val != 'nan' and val != 'NA': result[-1].append(float(val))
                    else: result[-1].append(np.nan)
                except:  # noqa: E722 (original)
                    result[-1].append(val)
    return result

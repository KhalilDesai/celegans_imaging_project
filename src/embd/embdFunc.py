"""
Port of development/embdFunc.py (Created on Sep 4, 2015, @author: Admin).

Changes: py3 syntax; bare ``raise`` in getNDDistance -> ValueError (DECISIONS R5); ``git`` is
imported lazily (only get_rnai_revision_version needs it); gene_key.csv path from fileLookup.
Arithmetic is unchanged.
"""
import numpy as np
import pandas as pd

from .fileLookup import FILENAME_GENE_NAMES, FILENAME_GENE_WB_INFO
from .myFunc import readCSV
from .myMath import getSig
from .varLookup import PARAM_NAMES

embsT = []


def getNDDistance(v1, v2, pNamesUse):
    ''' calculates distance between two vectors in N dimensional space. The vector component that are nan treated as dimensionality reduced '''

    if isinstance(v2, str) and v2 == 'control':
        v2 = np.ones_like(v1) * np.nan  # if v2 is control, set it to the origin
        for i in range(len(PARAM_NAMES)):  # iterates through the index of parameter names to analyze
            if PARAM_NAMES[i] in pNamesUse[0]:
                v2[i] = 0.
            if PARAM_NAMES[i] in pNamesUse[1]:
                v2[len(PARAM_NAMES) + i] = 0.
    if v1.size != v2.size:
        print('embdFunc.getNDDistance: wrong dimensionality of the supplied vectors')
        raise ValueError('embdFunc.getNDDistance: wrong dimensionality of the supplied vectors')
    N = len(pNamesUse[0]) + len(pNamesUse[1])
    vOverlap = np.invert(np.isnan(
        v1 * v2))  # overlapping vector, checks whether a parameter (coordinate) is not nan for each considered group (RNAi)
    if sum(vOverlap) == 0:
        return np.inf  # if there are no overlapping coordinates, returns infinity
    i1 = np.sum(np.invert(np.isnan(v1))) - np.sum(
        vOverlap)  # sums all coordinates that are not equal to nan for v1 and subtracts overlap
    i2 = np.sum(np.invert(np.isnan(v2))) - np.sum(
        vOverlap)  # sums all coordinates that are not equal to nan for v2 and subtracts overlap
    dist = np.sqrt(1. * N / (N - i1 - i2)) * np.linalg.norm((v1 - v2)[
                                                                vOverlap])  # imposes weighting to penalize for non-overlapping parameters and caluculates Euclidean Distance between 2 points
    return dist


def getNDAvgVector(vectors):
    ''' calculates average between list of vectors in N dimensional space. The vector component that are zero treated as dimensionality reduced '''
    vectors = np.array(vectors)  # gets list of numbers (vectors) from RNAi class and converts it to an array
    with np.errstate(invalid='ignore'):
        import warnings
        with warnings.catch_warnings():  # "Mean of empty slice" for all-NaN columns; result NaN as in numpy 1.13
            warnings.simplefilter('ignore', RuntimeWarning)
            avg = np.nanmean(vectors, axis=0)  # average of everything that is not nan
    avg[np.where(np.isnan(avg))] = np.nan  # fix...
    return avg


def get_interacting_genes(embdN, printflag=True):
    '''
    Retrieves a list of genes and evidence for interaction from CSV file- this information was mined from SimpleMine on wormbase 6/28/18
    '''
    gene = getGeneName(embdN)
    data = pd.read_csv(FILENAME_GENE_WB_INFO)
    gene_data = data[data.Sequence_Name == gene]
    interacting_genes = gene_data['Interacting_Gene'].tolist()
    ig_split = interacting_genes[0].split(',')
    if printflag: print('EMBD{0}: {1}, interacts with {2}'.format(embdN, getGeneCommonName(embdN), ig_split))
    return ig_split


def check_for_interaction(test_string, embdN):
    test = test_string
    interacting = get_interacting_genes(embdN, printflag=False)
    print('all interactions for {0}:\n {1}'.format(embdN, interacting))
    print('interactions matching string:')
    print([s for s in interacting if test in s])


def getGeneName(embdN):
    if isinstance(embdN, str): embdN = int(embdN[-4:])
    if embdN == 0: return 'control'
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == 'EMB_P{n:04}'.format(n=embdN))[0]
    return data[ind[0]][1]


def getGeneCommonName(embdN):
    if isinstance(embdN, str): embdN = int(embdN[-4:])
    if embdN == 0: return 'control'
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == 'EMB_P{n:04}'.format(n=embdN))[0]
    if ind.size > 0 and bool(data[ind[0]][2]):  # checks if anything is present in the cell
        return data[ind[0]][2]
    else:
        return 'unnamed'


def getEMBDNumber(geneName):
    if geneName == 'control': return 0
    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]
    dataA = np.array(data)
    ind = np.where(dataA == geneName)[0]
    if ind.size > 1:
        ind = ind[0]
    elif ind.size == 1:  # py2: `ind != [None]` (truthy for a 1-element array)
        return int(dataA[ind][0][0][-4:])
    else:
        return None


def get_param_names():
    """
    generates names for all parameters from varlookup.PARAM_NAMES
    :return: list of names
    """

    p_names = []
    for name in PARAM_NAMES:
        p_names.append(name + '_GLS')
        p_names.append(name + '_MS')
    return p_names


def split_params(params_all):
    """
    splits all params into lists of PARAM_NAMES GLS and PARAM_NAMES MS and the rest
    :param dict: dictionary of the parameters
    :return: three dictionaries: PARAM_NAMES GLS and PARAM_NAMES MS and the rest
    """
    params_gls = dict(zip(PARAM_NAMES, [np.nan for p in PARAM_NAMES]))
    params_ms = dict(zip(PARAM_NAMES, [np.nan for p in PARAM_NAMES]))
    other = {}
    for key in params_all:
        if key[-3:] == 'GLS' and key[:-4] in PARAM_NAMES:
            if params_all[key] is None:
                params_gls[key[:-4]] = np.nan
            else:
                params_gls[key[:-4]] = params_all[key]
        elif key[-3:] == '_MS' and key[:-3] in PARAM_NAMES:
            if params_all[key] is None:
                params_ms[key[:-3]] = np.nan
            else:
                params_ms[key[:-3]] = params_all[key]
        else:
            other[key] = params_all[key]
    return params_gls, params_ms, other


def get_rnai_revision_version():
    import git
    try:
        repo = git.Repo()
    except Exception:
        path = "../"
        repo = git.Repo(path=path)
    version = list(repo.iter_commits())[0]
    return str(version)[:7]


def make_dict_float(mydict):
    for key in mydict:
        if not isinstance(mydict[key], int):
            mydict[key] = float(np.round(mydict[key], getSig(mydict[key])+5))
    return mydict

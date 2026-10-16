"""
Port of development/PADassist.py: PAD-based ranking, PAD_RANK, best-match outputs.

Changes (docs/DECISIONS_AND_ISSUES.md, section P2):
* py3 syntax; 'Z:/...' paths -> FOLDER_IN; output file names are parameters whose default is the
  original date-stamped name;
* OUTPUT-NEUTRAL SPEED-UPS (maintainer request 2026-10-16): ``readCSV`` results of the four
  phenoDist files are memoised (the original re-read four 253k-row files on each of 503
  iterations), and ``PADFilter`` reuses the ``np.array(data)`` conversion of the same list
  object instead of rebuilding it on every call. Both are pure caches of unchanged inputs;
  results are identical (tests/test_port_regression.py);
* movies/images (show_best_match_movies, make_movies, convert_tifstack2mp4), Yanai-cluster
  lookups (CSV not public) and functions needing the missing dataValidation module are not ported.
"""
import time

import numpy as np

from ._compat import not_ported, py2_str
from .db_utils_embryos import get_columns, initialize_db
from .embdFunc import getEMBDNumber, getGeneCommonName
from .fileLookup import FILENAME_PAD
from .myFunc import readCSV as _readCSV_uncached
from .RNAiClass import RNAiClass
from .varLookup import FOLDER_IN

show_best_match_movies = not_ported('PADassist.show_best_match_movies')
make_movies = not_ported('PADassist.make_movies')
convert_tifstack2mp4 = not_ported('PADassist.convert_tifstack2mp4')
lookup_yanai_cluster = not_ported('PADassist.lookup_yanai_cluster')
check_yanai_cluster_by_EMBDgroup = not_ported('PADassist.check_yanai_cluster_by_EMBDgroup')
get_germlayer_yanai = not_ported('PADassist.get_germlayer_yanai')
findCliques = not_ported('PADassist.findCliques (needs the missing dataValidation module)')
getInteractions = not_ported('PADassist.getInteractions (needs IntegratedInteractome.csv)')

# --- output-neutral caches ---------------------------------------------------------------
_CSV_CACHE = {}
_ARRAY_CACHE = {}


def readCSV(fileName, delimiter=' '):
    """Memoised myFunc.readCSV (pure function of an unchanging file; callers only read the result)."""
    import os
    key = (fileName, delimiter, os.path.getmtime(fileName), os.path.getsize(fileName))
    if key not in _CSV_CACHE:
        _CSV_CACHE[key] = _readCSV_uncached(fileName, delimiter)
    return _CSV_CACHE[key]


def _as_array(data):
    """np.array(data) for the same (unmodified) list object is computed once."""
    k = id(data)
    hit = _ARRAY_CACHE.get(k)
    if hit is None or hit[0] is not data or hit[1] != len(data):
        hit = (data, len(data), np.array(data))
        _ARRAY_CACHE[k] = hit
    return hit[2]


def clear_caches():
    _CSV_CACHE.clear()
    _ARRAY_CACHE.clear()


# --- ranking ------------------------------------------------------------------------------
def PADFilter(rnaName, data, nTop):
    """
    Returns top ranked hits by pad
    :param rnaName: seed (EMBD string)
    :param data: data file
    :return: a list of lists for nTop hits. For a given seed gene (rnaName) each list consists of top hit embd name,
    PAD value and the number of overlapping dimensions('EMBD0076', 0.934496703417, 23.0)
    """
    if isinstance(rnaName, int):
        rnaName = 'EMBD{r:04}'.format(r=rnaName)
    nTop = nTop
    dataA = _as_array(data)  # original: np.array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[6] for p in dRNA]  # this is PAD
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][6], dRNA[i][4]) for i in inds[1:nTop + 1]]
    return res


def CSIFilter(rnaName, data):
    dataA = _as_array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[5] for p in dRNA]  # this is CSI
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][5], dRNA[i][4]) for i in inds[1:11]]
    return res


def distFilter(rnaName, data):
    dataA = _as_array(data)
    inds = np.where(dataA[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[2] for p in dRNA]  # this is Euclidean distance
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][2], dRNA[i][4]) for i in inds[1:11]]
    return res


def getPADFromFile(rn1, rn2, data):
    ind = np.where(data == 'EMBD{n:04}'.format(n=rn1))[0]
    if ind.size > 0:
        ind2 = np.where(data[ind] == 'EMBD{n:04}'.format(n=rn2))[0]
    else:
        print('!!!!!!!!!!!!!!!!!!!ALARM', rn1)
    if ind2.size > 1: ind2 = ind2[0]
    val = float(data[ind][ind2][6])
    return val


def getCSIFromFile(rn1, rn2, data):
    ind = np.where(data == 'EMBD{n:04}'.format(n=rn1))[0]
    ind2 = np.where(data[ind] == 'EMBD{n:04}'.format(n=rn2))[0]
    if ind2.size > 1: ind2 = ind2[0]
    val = float(data[ind][ind2][5])
    return val


def makeNetwork(edgeData):
    '''loads genes as nodes and PAD values as edges using networkx'''
    import networkx as nx
    G = nx.Graph()
    for i in range(1, 504):
        G.add_node(i)
    th = 0.74
    edgeData = np.array(edgeData)
    inds = np.where(edgeData[:, 6].astype(float) > th)
    edgeData = edgeData[inds]
    for row in edgeData:
        n1 = int(row[0][-4:])  # node1
        n2 = int(row[1][-4:])  # node2
        weight = float(row[6])  # PAD
        if n1 != n2:
            G.add_edge(n1, n2, weight=weight)  # adds edges
    return G  # graph of all genes with PAD as edge value


def getPairsParams(seed, data):
    '''get PAD values and most correlated parameter names between seed and top ten neighbors (i.e. 156)'''

    rnaSeed = RNAiClass(seed)
    rnaSeed.setNDPosition()
    pairList = PADFilter('EMBD{r:04}'.format(r=seed), data, 10)
    rList = [RNAiClass(int(p[0][5:])) for p in pairList]
    print('Pairs for EMBD{r:04}'.format(r=seed))
    i = 0
    for p in pairList:
        print('{r}: PAD={pad}, nOverlap={o}, {par}'.format(r=p[0], pad=rnaSeed.getPAD(rList[i]), o=p[2],
                                                           par=rnaSeed.getSimParams(rList[i])))
        i += 1


def getCommonPairs(seedList, data):
    allPairList = []
    allPADvals = []
    for seed in seedList:
        topGenes = PADFilter('EMBD{r:04}'.format(r=seed), data, 10)
        topGenes = np.array(topGenes)
        allPairList.append(topGenes[:, 0])
        allPADvals.append(topGenes[:, 1])
    uniqueList, countList = np.unique(allPairList, return_counts=True)
    allPairList = np.array(allPairList)
    allPADvals = np.array(allPADvals)
    uniquePADvals = []
    for gene in uniqueList:
        uniquePADvals.append(np.mean(allPADvals[np.where(allPairList == gene)].astype(float)))
    index = np.lexsort((uniquePADvals, countList))
    for i in index[::-1]:
        print(uniqueList[i], countList[i], uniquePADvals[i])
    return [int(uniqueList[i][-4:]) for i in index[::-1]]


def get_best_match(RNAi, data):
    '''input- RNAi object- i.e EMBD0076
    Finds top matches using PAD (orCSI) by calling a CSV file, returns a list of ints [12,33,126] representing the genes with the closest phenotypes'''

    seed = RNAi[-4:]
    seedObj = RNAiClass(int(seed))  # noqa: F841 (original)
    top_hits = [seed]
    top_hit_data = np.array(PADFilter(RNAi, data, 5))[:, 0]
    top_hit_vals = np.array(PADFilter(RNAi, data, 5))[:, 1]
    for top_hit in top_hit_data:
        top_hits.append(top_hit[-4:])
    print(top_hits, seed, top_hit_vals)
    return top_hits, seed, top_hit_vals


# --- data access --------------------------------------------------------------------------
def get_data():
    """
    Retrieves all pairwise information for sets of genes (from csv; the mysql variant is kept below).
    """
    return get_data_from_csv()


def get_data_from_mysql():
    conn, cursor = initialize_db()
    sql = "SELECT 'EMBD' || substr('0000' || rnai1, -4, 4), 'EMBD' || substr('0000' || rnai2, -4, 4), dist, " \
          "dist2wt, dims, CSI, PAD FROM pair_lookup"  # original MySQL: CONCAT('EMBD', LPAD(rnai1, 4, '0')), ...
    cursor.execute(sql)
    dataset = list(cursor.fetchall())
    return dataset


def get_data_from_csv():
    dataset = readCSV(FILENAME_PAD, ',')[3:]
    return dataset


def update_pad_rank(cursor, conn, data, embd_ids=range(1, 504)):
    """
    adds genes to the pad_rank database
    """
    from .db_utils_embryos import insert_row, update_row
    table = 'pad_rank'
    columns = get_columns(cursor, table)  # retrieves columns from table as list
    for i in embd_ids:
        print(i)
        top_hits = PADFilter(i, data, 51)
        count = 0
        colms = {}
        for col in columns[1:]:
            colms[col] = 0  # add a default value to populate key value pair
        for hit in top_hits:
            rnai = i
            hit_name = int(hit[0][4:])
            pad = hit[1]
            count += 1
            rank_col = 'rank'+'{}'.format(count)
            pad_col = 'pad'+'{}'.format(count)
            colms['rnai'] = rnai
            colms[rank_col] = hit_name
            colms[pad_col] = pad

        sql = "SELECT id FROM {t} WHERE rnai={r}".format(t=table, r=i)
        cursor.execute(sql)
        vals = cursor.fetchall()
        if len(vals) > 0:
            update_row(int(vals[0][0]), colms, cursor, table)
        elif len(vals) == 0:
            insert_row(colms, cursor, table)
    conn.commit()
    conn.close()


# --- best-match outputs -------------------------------------------------------------------
PHENODIST_FILES = {'opt': 'phenoDist_rank_04142022.csv', 'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
                   'ms': 'phenoDist_rank_MS_ONLY05102022.csv', 'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}


def _load_all_sets():
    return (readCSV(FOLDER_IN + PHENODIST_FILES['opt'], ',')[3:],
            readCSV(FOLDER_IN + PHENODIST_FILES['gls'], ',')[3:],
            readCSV(FOLDER_IN + PHENODIST_FILES['ms'], ',')[3:],
            readCSV(FOLDER_IN + PHENODIST_FILES['all'], ',')[3:])


def output_all_best_matches_to_csv_combined(fileName=None, embd_ids=range(1, 504)):
    '''
    for all conditions, finds top matches using PAD, writes to a CSV file: shows top matches with alternate parameter sets'''
    if fileName is None:
        fileName = FOLDER_IN + 'all_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y'))
    fd = open(fileName, 'w', newline='')
    for i in embd_ids:
        RNAi = 'EMBD{r:04}'.format(r=i)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))  # noqa: F841 (original)

        data_opt, data_GLS_params_only, data_MS_params_only, data_all_params = _load_all_sets()

        top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]
        top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]
        top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]
        top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]
        top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]
        top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]
        top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]
        top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]
        fd.write('rank, query EMBD_ID, query genename, hit EMBD_ID - optimal parameters, hit genename - optimal parameters, hit PAD value - optimal parameters, hit EMBD_ID - GLS parameters, hit genename - GLS parameters, hit PAD value - GLS parameters, hit EMBD_ID - MS parameters, hit genename - MS parameters, hit PAD value - MS parameters+ \n')
        all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
        for i in range(len(top_hit_data_opt)):  # noqa: B020 (shadows the outer i, as in the original)
            rank = i + 1
            query = seed
            query_name = getGeneCommonName(seed)
            opt_hit = top_hit_data_opt[i][-4:]
            opt_hit_name = getGeneCommonName(opt_hit)
            opt_hit_val = top_hit_vals_opt[i]
            GLS_hit = top_hit_data_GLS[i][-4:]
            GLS_hit_name = getGeneCommonName(GLS_hit)
            GLS_hit_val = top_hit_vals_GLS[i]
            MS_hit = top_hit_data_MS[i][-4:]
            MS_hit_name = getGeneCommonName(MS_hit)
            MS_hit_val = top_hit_vals_MS[i]
            all_hit = top_hit_data_all_params[i][-4:]
            all_hit_name = getGeneCommonName(all_hit)
            all_hit_val = top_hit_vals_all_params[i]
            all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val, all_hit, all_hit_name, all_hit_val])
            fd.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val))
            if i == 9:
                fd.write('\n')
    fd.close()
    return fileName


def output_all_best_matches_to_one_csv_combined(fileName=None, embd_ids=range(1, 504)):
    '''
    for all conditions, finds top matches using PAD, writes to one CSV file for all conditions: shows top matches and PAD (hit-1: 0.765) with alternate parameter sets'''
    import pandas as pd
    all_hits_matrix = [['query EMBD_ID', 'query genename', 'hits - optimal parameters', 'hits - GLS only parameters', 'hits - MS only parameters',
                        'hits - all parameters']]

    for i in embd_ids:
        RNAi = 'EMBD{r:04}'.format(r=i)
        print(RNAi)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))  # noqa: F841 (original)
        query = seed
        query_name = getGeneCommonName(seed)

        data_opt, data_GLS_params_only, data_MS_params_only, data_all_params = _load_all_sets()

        sets = []  # the original repeats this block four times verbatim; identical logic
        for data in (data_opt, data_GLS_params_only, data_MS_params_only, data_all_params):
            top_hit_data = np.array(PADFilter(RNAi, data, 10))[:, 0]
            top_hit_data_gene_names = [getGeneCommonName(top_hit_data[j][-4:]) for j in range(len(top_hit_data))]
            top_hit_vals = np.array(PADFilter(RNAi, data, 10))[:, 1]
            top_hit_set = []
            for k in range(10):
                val = float(top_hit_vals[k])
                if val > 0.50:
                    val = py2_str(val)[0:5]  # py2 str(float); DECISIONS B6
                    top_hit_set.append(top_hit_data_gene_names[k] + ': ' + val)
            sets.append(top_hit_set)

        all_hits_matrix.append([query, query_name] + sets)

    if fileName is None:
        fileName = FOLDER_IN + 'hits_matrix_v2_{d}.csv'.format(d=time.strftime('%m%d%Y'))
    pd.DataFrame(all_hits_matrix).to_csv(fileName)
    return fileName


def output_best_matches_to_csv(RNAi):
    raise NotImplementedError('PADassist.output_best_matches_to_csv: per-gene variant writing into the '
                              'web-content tree; not on the published-output path, not ported')

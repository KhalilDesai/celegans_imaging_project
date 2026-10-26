"""
Port of development/PADassist.py: PAD-based ranking, PAD_RANK, best-match outputs.

Changes (docs/DECISIONS_AND_ISSUES.md, section P2):
* py3 syntax; 'Z:/...' paths -> FOLDER_IN; output file names are parameters whose default is the
  original date-stamped name;
* OUTPUT-NEUTRAL SPEED-UPS (maintainer request 2026-10-16): ``readCSV`` results of the four
  phenoDist files are memoised (the original re-read four 253k-row files on each of 503
  iterations), and ``PADFilter``/``CSIFilter``/``distFilter`` find a seed's rows through a label index built
  once per list instead of converting the whole list to a string array on every call. Both are pure caches of unchanged inputs;
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


def _rows_for(data, rnaName):
    """Indices of the rows of ``data`` whose first field equals ``rnaName``, ascending.

    Identical to the original ``np.where(np.array(data)[:, 0] == rnaName)[0]`` but uses a label->rows
    index built once per (unchanged) list object instead of converting the whole 253k x 7 list to a
    string array on every call (that conversion was ~230 MB per file in numpy 2 and caused an
    out-of-memory kill on an 8 GB machine).
    """
    k = id(data)
    hit = _ARRAY_CACHE.get(k)
    if hit is None or hit[0] is not data or hit[1] != len(data):
        index = {}
        for i, row in enumerate(data):
            index.setdefault(str(row[0]), []).append(i)
        hit = (data, len(data), {lab: np.array(ix, dtype=np.intp) for lab, ix in index.items()})
        _ARRAY_CACHE[k] = hit
    return hit[2].get(rnaName, np.array([], dtype=np.intp))


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
    inds = _rows_for(data, rnaName)  # original: np.where(np.array(data)[:, 0] == rnaName)[0]
    dRNA = [data[i] for i in inds]
    pads = [p[6] for p in dRNA]  # this is PAD
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][6], dRNA[i][4]) for i in inds[1:nTop + 1]]
    return res


def CSIFilter(rnaName, data):
    inds = _rows_for(data, rnaName)
    dRNA = [data[i] for i in inds]
    pads = [p[5] for p in dRNA]  # this is CSI
    inds = np.argsort(pads)[::-1]
    res = [(dRNA[i][1], dRNA[i][5], dRNA[i][4]) for i in inds[1:11]]
    return res


def distFilter(rnaName, data):
    inds = _rows_for(data, rnaName)
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


#     for d in data:
#         if d[0]=='EMBD{n:04}'.format(n=rn1) and d[1]=='EMBD{n:04}'.format(n=rn2): return d[6]
#         elif d[0]=='EMBD{n:04}'.format(n=rn2) and d[1]=='EMBD{n:04}'.format(n=rn1): return d[6]
#     return np.nan


def makeNetwork(edgeData):
    '''loads genes as nodes and PAD values as edges using networkx'''
    import networkx as nx  # PORT: was a module-level import
    G = nx.Graph()  # creates empty graph with no nodes or edges, can assign graph attributes i.e. G = nx.Graph(day="Friday")
    for i in range(1, 504):
        G.add_node(i)
    #     G.add_nodes_from([range(1,504)], color = 'blue') #add a list of nodes
    th = 0.74
    #     th2= 0.7
    edgeData = np.array(edgeData)
    inds = np.where(edgeData[:, 6].astype(float) > th)
    #     inds = getSubIndex(edgeData[:,6].astype(float), th2, th)
    edgeData = edgeData[inds]
    for row in edgeData:
        n1 = int(row[0][-4:])  # node1
        n2 = int(row[1][-4:])  # node2
        weight = float(row[6])  # PAD
        #         weight = float(row[5])#CSI
        if n1 != n2:
            G.add_edge(n1, n2, weight=weight)  # adds edges
        #     for i in range(1,3):
        #         print G[i]
    return G  # graph of all genes with PAD as edge value


def getPairsParams(seed, data):
    '''get PAD values and most correlated parameter names between seed and top ten neighbors (i.e. 156)'''

    rnaSeed = RNAiClass(seed)
    rnaSeed.setNDPosition()
    pairList = PADFilter('EMBD{r:04}'.format(r=seed),
                         data, 10)  # returns a list of top hits by PAD values- each as gene, PAD, dimensions
    rList = [RNAiClass(int(p[0][5:])) for p in pairList]
    print('Pairs for EMBD{r:04}'.format(r=seed))
    i = 0
    for p in pairList:
        print('{r}: PAD={pad}, nOverlap={o}, {par}'.format(r=p[0], pad=rnaSeed.getPAD(rList[i]), o=p[2],
                                                           par=rnaSeed.getSimParams(rList[i])))
        i += 1


def getCommonPairs(seedList, data):
    allPairList = []  # list of EMBD numbers that are paired with seed
    allPADvals = []  # list of PAD values that from genes paired with seed
    for seed in seedList:
        topGenes = PADFilter('EMBD{r:04}'.format(r=seed),
                             data, 10)  # returns list of lists- each list has gene, PAD and number of dimensions
        topGenes = np.array(topGenes)
        allPairList.append(topGenes[:, 0])  # appends genes
        allPADvals.append(topGenes[:, 1])  # appends PAD vals
    uniqueList, countList = np.unique(allPairList,
                                      return_counts=True)  # returns two lists- first list is all unique genes, second list is count for how many times genes are in list
    allPairList = np.array(allPairList)
    allPADvals = np.array(allPADvals)
    uniquePADvals = []  # list of PAD values for list of genes (reduced to unique list)
    for gene in uniqueList:
        uniquePADvals.append(np.mean(allPADvals[np.where(allPairList == gene)].astype(
            float)))  # finds the TRUE indicies where gene is in allPairList, retrieves the PAD values from those indices. If there are more than one pad value, appends the average
    index = np.lexsort((uniquePADvals,
                        countList))  # returns a list of indices corresponded to sorting by countList, then by uniquePADvals
    for i in index[::-1]:  # [::-1] slices through the list in reverse order
        print(uniqueList[i], countList[i], uniquePADvals[i])
    return [int(uniqueList[i][-4:]) for i in index[::-1]]


#     np.argsort(uniqueList[:,1])


def get_best_match(RNAi, data):
    '''input- RNAi object- i.e EMBD0076
    Finds top matches using PAD (orCSI) by calling a CSV file, returns a list of ints [12,33,126] representing the genes with the closest phenotypes'''

    seed = RNAi[-4:]
    seedObj = RNAiClass(int(seed))
    top_hits = [
        seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data = np.array(PADFilter(RNAi, data, 5))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals = np.array(PADFilter(RNAi, data, 5))[:, 1]  # retrieves the PAD values for top matches

    ## top_hit_data = np.array(CSIFilter(RNAi, data))[:,0] #retrieves the top matches i.e. 'EMBD0076' based on CSI
    ## top_hit_vals = np.array(CSIFilter(RNAi, data))[:,1] #retrieves the CSI values for top matches

    for top_hit in top_hit_data:
        top_hits.append(top_hit[-4:])  # slices the last four digits of string to remove the 'EMBD'
    print(top_hits, seed, top_hit_vals)
    return top_hits, seed, top_hit_vals
def get_data():
    """
    Retrieves all pairwise information for sets of genes. for each gene pair: (rnai1 label, rnai2 label, distance,
    distance to wt, overlapping dimensions, csi and pad from either mysql or csv

    :return: data: a list of lists. Each lists consists of rnai1 label, rnai2 label, distance,
    distance to wt, overlapping dimensions, csi and pad
    """

    '''get data from mysql'''
    # return get_data_from_mysql()

    # '''get data from csv'''
    return get_data_from_csv()


def get_data_from_mysql():
    conn, cursor = initialize_db()
    # PORT: the SQLite emulation has no CONCAT/LPAD; original: SELECT CONCAT('EMBD', LPAD(rnai1, 4, '0')), ...
    sql = "SELECT 'EMBD' || substr('0000' || rnai1, -4, 4), 'EMBD' || substr('0000' || rnai2, -4, 4), dist, " \
          "dist2wt, dims, CSI, PAD FROM pair_lookup"  # LPAD Left PADs the INT rnai1 value with 'EMBD' and zeros, up to 4 digits

    cursor.execute(sql)
    dataset = list(cursor.fetchall())
    return dataset


def get_data_from_csv():
    dataset = readCSV(FILENAME_PAD, ',')[3:]
    return dataset


def update_pad_rank(cursor, conn, data, embd_ids=range(1, 504)):
    """
    adds genes to the pad_rank database
    :param cursor: mysqldb cursor
    :data: data set containing gene pairs, PAD, distance, CSI values generated by get_data(), which
    either populates data from csv or mysql
    :return: None
    """
    from .db_utils_embryos import update_row, insert_row
    table = 'pad_rank'
    columns = get_columns(cursor, table)  # retrieves columns from table as list
    for i in embd_ids:  # PORT: the original range(1, 504) is the default
        print(i)
        top_hits = PADFilter(i, data, 51)
        count = 0
        colms = {}
        for col in columns[1:]:
            colms[col] = 0 #add a default value to populate key value pair
        for hit in top_hits:
            rnai = i
            hit_name = int(hit[0][4:])
            pad = hit[1]
            count+=1
            rank_col = 'rank'+'{}'.format(count)
            pad_col = 'pad'+'{}'.format(count)
            colms['rnai']=rnai
            colms[rank_col]=hit_name
            colms[pad_col]=pad

        sql = "SELECT id FROM {t} WHERE rnai={r}".format(t=table, r=i)
        cursor.execute(sql)
        vals = cursor.fetchall()
        if len(vals) > 0:
            update_row(int(vals[0][0]), colms, cursor, table)
        elif len(vals) == 0:
            insert_row(colms, cursor, table)
    conn.commit()
    conn.close()
PHENODIST_FILES = {'opt': 'phenoDist_rank_04142022.csv', 'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
                   'ms': 'phenoDist_rank_MS_ONLY05102022.csv', 'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}


def output_all_best_matches_to_csv_combined(fileName=None, embd_ids=range(1, 504)):

    '''
    for all conditions, finds top matches using PAD, writes to a CSV file: shows top matches with alternate parameter sets'''
    import time

    if fileName is None:  # PORT: was 'Z:\\all_best_matches_<date>.csv'
        fileName = FOLDER_IN + 'all_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y'))
    fd = open(fileName, 'w', newline='')
    for i in embd_ids:  # PORT: the original range(1,504) is the default
        RNAi = 'EMBD{r:04}'.format(r=i)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))

        data_opt = readCSV(FOLDER_IN + 'phenoDist_rank_04142022.csv', ',')[3:]
        data_GLS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
        data_MS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
        data_all_params = readCSV(FOLDER_IN + 'phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


        # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_GLS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_MS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

        # top_hits_all_params = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches
        fd.write('rank, query EMBD_ID, query genename, hit EMBD_ID - optimal parameters, hit genename - optimal parameters, hit PAD value - optimal parameters, hit EMBD_ID - GLS parameters, hit genename - GLS parameters, hit PAD value - GLS parameters, hit EMBD_ID - MS parameters, hit genename - MS parameters, hit PAD value - MS parameters+ \n')
        all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
        for i in range(len(top_hit_data_opt)):
            rank = i + 1
            query = seed
            query_name = getGeneCommonName(seed)
            opt_hit =top_hit_data_opt[i][-4:]
            opt_hit_name = getGeneCommonName(opt_hit)
            opt_hit_val = top_hit_vals_opt[i]
            GLS_hit =top_hit_data_GLS[i][-4:]
            GLS_hit_name = getGeneCommonName(GLS_hit)
            GLS_hit_val = top_hit_vals_GLS[i]
            MS_hit =top_hit_data_MS[i][-4:]
            MS_hit_name = getGeneCommonName(MS_hit)
            MS_hit_val = top_hit_vals_MS[i]
            all_hit =top_hit_data_all_params[i][-4:]
            all_hit_name = getGeneCommonName(all_hit)
            all_hit_val = top_hit_vals_all_params[i]
            all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])
            fd.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val))
            if i==9:
                fd.write('\n')

            # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
            # hit_name = getGeneCommonName(top_hit_opt[-4:])
        # all_hits_matrix = np.array(all_hits_matrix)
        # print(all_hits_matrix)


        # import pandas as pd
        # savepath = 'Z:/WebContent/{0}/{0}_hits_matrix.csv'.format(RNAi)
        # pd.DataFrame(all_hits_matrix).to_csv(savepath)
    fd.close()
    return fileName  # PORT
        # return(all_hits_matrix)
        # print top_hits, seed, top_hit_vals
        # return top_hits, seed, top_hit_vals


def output_all_best_matches_to_one_csv_combined(fileName=None, embd_ids=range(1, 504)):
    '''
    IN PROGRESS__for all conditions, finds top matches using PAD, writes to one CSV file for all conditions: shows top matches  and PAD (hit-1: 0.765) with alternate parameter sets'''
    import time

    # all_hits_matrix = [['query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters',
    #                     'hit genename - optimal parameters', 'hit PAD value - optimal parameters',
    #                     'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters',
    #                     'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters',
    #                     'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters',
    #                     'hit genename - all parameters', 'hit PAD value - all parameters']]

    all_hits_matrix = [['query EMBD_ID', 'query genename', 'hits - optimal parameters', 'hits - GLS only parameters', 'hits - MS only parameters',
                        'hits - all parameters']]

    # fd = open('Z:\\Automated_analysis\\combined_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y')),
    #           'w')
    # fd.write(
    #     'query EMBD_ID, query genename, hit EMBD_ID - optimal parameters, hit genename - optimal parameters, hit PAD value - optimal parameters, hit EMBD_ID - GLS parameters, hit genename - GLS parameters, hit PAD value - GLS parameters, hit EMBD_ID - MS parameters, hit genename - MS parameters, hit PAD value - MS parameters, hit EMBD_ID - All parameters, hit genename - All parameters, hit PAD value - All parameters+ \n')

    for i in embd_ids:  # PORT: the original range(1,504) is the default
        RNAi = 'EMBD{r:04}'.format(r=i)
        print(RNAi)
        seed = RNAi[-4:]
        seedObj = RNAiClass(int(seed))
        query = seed
        query_name = getGeneCommonName(seed)

        data_opt = readCSV(FOLDER_IN + 'phenoDist_rank_04142022.csv', ',')[3:]
        data_GLS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
        data_MS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
        data_all_params = readCSV(FOLDER_IN + 'phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


        # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
        top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        # top_hit_data_opt_gene_names = []

        # for j in range(len(top_hit_data_opt)):
        #     hit_name = getGeneCommonName(top_hit_data_opt[j][-4:])
        #     top_hit_data_opt_gene_names.append(hit_name)

        top_hit_data_opt_gene_names= [getGeneCommonName(top_hit_data_opt[j][-4:]) for j in range(len(top_hit_data_opt))]
        top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_opt_set = []
        for k in range(10):
            val = float(top_hit_vals_opt[k])
            if val>0.50:
                val = py2_str(val)[0:5]  # PORT: py2 str(float) (B6)
                top_hit_opt_set.append(top_hit_data_opt_gene_names[k]+': '+ val)

        top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_GLS_gene_names = [getGeneCommonName(top_hit_data_GLS[j][-4:]) for j in
                                       range(len(top_hit_data_GLS))]
        top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_GLS_set = []
        for l in range(10):
            val = float(top_hit_vals_GLS[l])
            if val > 0.50:
                val = py2_str(val)[0:5]  # PORT: py2 str(float) (B6)
                top_hit_GLS_set.append(top_hit_data_GLS_gene_names[l] + ': ' + val)

        top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_MS_gene_names = [getGeneCommonName(top_hit_data_MS[j][-4:]) for j in
                                       range(len(top_hit_data_MS))]
        top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_MS_set = []
        for m in range(10):
            val = float(top_hit_vals_MS[m])
            if val > 0.50:
                val = py2_str(val)[0:5]  # PORT: py2 str(float) (B6)
                top_hit_MS_set.append(top_hit_data_MS_gene_names[m] + ': ' + val)

        top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
        top_hit_data_all_gene_names = [getGeneCommonName(top_hit_data_all_params[j][-4:]) for j in
                                       range(len(top_hit_data_all_params))]
        top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches
        top_hit_all_set = []
        for n in range(10):
            val = float(top_hit_vals_all_params[n])
            if val > 0.50:
                val = py2_str(val)[0:5]  # PORT: py2 str(float) (B6)
                top_hit_all_set.append(top_hit_data_all_gene_names[n] + ': ' + val)

        # fd.write(
        #     '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(query, query_name, top_hit_data_opt.tolist(), top_hit_data_opt_gene_names,
        #                                                                  top_hit_vals_opt.tolist(), top_hit_data_GLS.tolist(), top_hit_data_GLS_gene_names,
        #                                                                  top_hit_vals_GLS.tolist(), top_hit_data_MS.tolist(), top_hit_data_MS_gene_names, top_hit_vals_MS.tolist(), top_hit_data_all_params.tolist(), top_hit_data_all_gene_names, top_hit_vals_all_params.tolist()))

        # all_hits_matrix.append(
        #     [query, query_name, top_hit_data_opt.tolist(),top_hit_data_opt_gene_names,top_hit_vals_opt.tolist(),top_hit_data_GLS.tolist(),top_hit_data_GLS_gene_names,top_hit_vals_GLS.tolist(),top_hit_data_MS.tolist(),top_hit_data_MS_gene_names, top_hit_vals_MS.tolist(),top_hit_data_all_params.tolist(),top_hit_data_all_gene_names,top_hit_vals_all_params.tolist()])

        all_hits_matrix.append([query, query_name, top_hit_opt_set, top_hit_GLS_set, top_hit_MS_set,top_hit_all_set])

    import pandas as pd
    savepath = fileName if fileName is not None else \
        FOLDER_IN + 'hits_matrix_v2_{d}.csv'.format(d=time.strftime('%m%d%Y'))  # PORT: was 'Z:/'
    pd.DataFrame(all_hits_matrix).to_csv(savepath)
    return savepath  # PORT



        # all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit
        # all_hits_matrix = [['rank', 'query EMBD_ID', 'query genename', 'hit EMBD_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit EMBD_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit EMBD_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit EMBD_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
        # for i in range(len(top_hit_data_opt)):
        #     rank = i + 1
        #     query = seed
        #     query_name = getGeneCommonName(seed)
        #     opt_hit =top_hit_data_opt[i][-4:]
        #     opt_hit_name = getGeneCommonName(opt_hit)
        #     opt_hit_val = top_hit_vals_opt[i]
        #     GLS_hit =top_hit_data_GLS[i][-4:]
        #     GLS_hit_name = getGeneCommonName(GLS_hit)
        #     GLS_hit_val = top_hit_vals_GLS[i]
        #     MS_hit =top_hit_data_MS[i][-4:]
        #     MS_hit_name = getGeneCommonName(MS_hit)
        #     MS_hit_val = top_hit_vals_MS[i]
        #     all_hit =top_hit_data_all_params[i][-4:]
        #     all_hit_name = getGeneCommonName(all_hit)
        #     all_hit_val = top_hit_vals_all_params[i]
        #     all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])
        #     fd.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n'.format(rank, query, query_name, opt_hit, opt_hit_name, opt_hit_val, GLS_hit, GLS_hit_name, GLS_hit_val, MS_hit, MS_hit_name, MS_hit_val))
        #     if i==9:
        #         fd.write('\n')

            # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
            # hit_name = getGeneCommonName(top_hit_opt[-4:])
        # all_hits_matrix = np.array(all_hits_matrix)
        # print(all_hits_matrix)


        # import pandas as pd
        # savepath = 'Z:/WebContent/{0}/{0}_hits_matrix.csv'.format(RNAi)
        # pd.DataFrame(all_hits_matrix).to_csv(savepath)
    # fd.close()
        # return(all_hits_matrix)
        # print top_hits, seed, top_hit_vals
        # return top_hits, seed, top_hit_vals


def output_best_matches_to_csv(RNAi):
    '''input- RNAi object- i.e EMBD0076
    Finds top matches using PAD (orCSI) by calling a CSV file, returns a list of ints [12,33,126] representing the genes with the closest phenotypes'''

    seed = RNAi[-4:]
    seedObj = RNAiClass(int(seed))

    data_opt = readCSV(FOLDER_IN + 'phenoDist_rank_04142022.csv', ',')[3:]
    data_GLS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_GLS_ONLY05102022.csv', ',')[3:]
    data_MS_params_only = readCSV(FOLDER_IN + 'phenoDist_rank_MS_ONLY05102022.csv', ',')[3:]
    data_all_params = readCSV(FOLDER_IN + 'phenoDist_rank_ALL_PARAMS05102022.csv', ',')[3:]


    # top_hits_opt = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_opt = np.array(PADFilter(RNAi, data_opt, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_GLS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_GLS = np.array(PADFilter(RNAi, data_GLS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_MS = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_MS = np.array(PADFilter(RNAi, data_MS_params_only, 10))[:, 1]  # retrieves the PAD values for top matches

    # top_hits_all_params = [seed]  # top hits list is seeded with the RNAi of interest and is populated with top matches (slices the last four digits of string to remove the 'EMBD')
    top_hit_data_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 0]  # retrieves the top matches i.e. 'EMBD0076'
    top_hit_vals_all_params = np.array(PADFilter(RNAi, data_all_params, 10))[:, 1]  # retrieves the PAD values for top matches

    all_hits_matrix = [['rank', 'query unique_ID', 'query genename', 'hit unique_ID - optimal parameters', 'hit genename - optimal parameters', 'hit PAD value - optimal parameters', 'hit unique_ID - GLS parameters', 'hit genename - GLS parameters', 'hit PAD value - GLS parameters', 'hit unique_ID - MS parameters', 'hit genename - MS parameters', 'hit PAD value - MS parameters', 'hit unique_ID - all parameters', 'hit genename - all parameters', 'hit PAD value - all parameters']]
    for i in range(len(top_hit_data_opt)):
        rank = i + 1
        query = seed
        query_name = getGeneCommonName(seed)
        opt_hit =top_hit_data_opt[i][-4:]
        opt_hit_name = getGeneCommonName(opt_hit)
        opt_hit_val = top_hit_vals_opt[i]
        GLS_hit =top_hit_data_GLS[i][-4:]
        GLS_hit_name = getGeneCommonName(GLS_hit)
        GLS_hit_val = top_hit_vals_GLS[i]
        MS_hit =top_hit_data_MS[i][-4:]
        MS_hit_name = getGeneCommonName(MS_hit)
        MS_hit_val = top_hit_vals_MS[i]
        all_hit =top_hit_data_all_params[i][-4:]
        all_hit_name = getGeneCommonName(all_hit)
        all_hit_val = top_hit_vals_all_params[i]
        all_hits_matrix.append([rank, query, query_name, opt_hit, opt_hit_name,opt_hit_val, GLS_hit, GLS_hit_name,GLS_hit_val, MS_hit, MS_hit_name,MS_hit_val, all_hit, all_hit_name,all_hit_val])



        # all_hits_matrix.append(top_hit_data_opt[i][-4:])  # slices the last four digits of string to remove the 'EMBD'
        # hit_name = getGeneCommonName(top_hit_opt[-4:])
    all_hits_matrix = np.array(all_hits_matrix)
    print(all_hits_matrix)

    import pandas as pd
    savepath = FOLDER_IN + '{0}_hits_matrix.csv'.format(RNAi)  # PORT: was 'Z://'
    pd.DataFrame(all_hits_matrix).to_csv(savepath)

    return(all_hits_matrix)
    # print top_hits, seed, top_hit_vals
    # return top_hits, seed, top_hit_vals

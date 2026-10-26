"""Rebuild hand-ported functions from the ORIGINAL source so that their comments are preserved.

For every (port module, original file, function, edits) entry below, the function's original
source (body, docstring and comments, including commented-out code) is copied and only the listed
edits are applied. The result replaces that function in the port module. Each edit anchor must
match the stated number of times (default once), so the edit list is a complete, auditable record
of every change made in the port.

Run: python tools/port_functions.py        (idempotent; then run the test suite)
"""
import os
import re

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORT = 'PORT'  # marker used in comments added by the port

DBG = 'development/db_utils_genes.py'
DBE = 'embryos_analyze/db_utils_embryos.py'
PDC = 'development/phenoDistCalc.py'
PAD = 'development/PADassist.py'
RNA = 'development/RNAiClass.py'
EMB = 'development/embdFunc.py'
MTH = 'embryos_analyze/myMath.py'

CONN = "    conn = _conn_of(cursor, conn)  # PORT: the original used the __main__-level global `conn`\n"


def zpath(name):
    return ("readCSV('Z:/%s', ',')" % name, "readCSV(FOLDER_IN + '%s', ',')" % name)


PHENODIST_READS = [zpath(n) for n in ('phenoDist_rank_04142022.csv', 'phenoDist_rank_GLS_ONLY05102022.csv',
                                      'phenoDist_rank_MS_ONLY05102022.csv', 'phenoDist_rank_ALL_PARAMS05102022.csv')]

SPEC = {
    'src/embd/db_utils_genes.py': [
        (DBG, 'get_columns', []),
        (DBG, 'get_tot_num_genes', []),
        (DBG, 'get_id_from_embd', [("isinstance(embd_id, basestring):", "isinstance(embd_id, str):  # PORT: py2 basestring")]),
        (DBG, 'get_id_from_genename', []),
        (DBG, 'get_max_id', []),
        (DBG, 'create_columns_pad_rank', []),
        (DBG, 'get_reciprocal_matches', []),
        (DBG, 'update_automated_description', [("    from groupAssesment import", "    from .groupAssesment import")]),
        (DBG, 'update_top_hits_PAD', [
            ("def update_top_hits_PAD(cursor):", "def update_top_hits_PAD(cursor, conn=None, embd_ids=ORIGINAL_TOP_HITS_LIST):"),
            ("    from embdFunc import getGeneCommonName\n", "    from .embdFunc import getGeneCommonName\n" + CONN),
            ("    for i in [1406,674,952,866,885,703,1360,1904,1903,1905]:\n",
             "    for i in embd_ids:  # PORT: the hand-edited list [1406,674,952,866,885,703,1360,1904,1903,1905] is the default\n"),
            ("        print('starting {0}').format(i)\n", "        print('starting {0}'.format(i))  # PORT: py2 print statement\n")]),
        (DBG, 'update_WT_like_or_phenotype', [
            ("def update_WT_like_or_phenotype(cursor):", "def update_WT_like_or_phenotype(cursor, conn=None):"),
            ("    for i in range(1,504):\n        print('starting {0}').format(i)\n",
             CONN + "    for i in range(1,504):\n        print('starting {0}'.format(i))  # PORT: py2 print statement\n")]),
        (DBG, 'update_column_reciprocal_matches', [
            ("def update_column_reciprocal_matches(cursor):",
             "def update_column_reciprocal_matches(cursor, conn=None, embd_ids=range(310, 311)):"),
            ("    from embdFunc import getGeneCommonName\n", "    from .embdFunc import getGeneCommonName\n" + CONN),
            ("    for i in range(310,311):\n", "    for i in embd_ids:  # PORT: the hand-edited range(310,311) is the default\n")]),
        (DBG, 'get_sql_data_for_man_group_conditions', [
            ("    df = pd.read_sql_query(sql, conn)",
             "    df = pd.DataFrame(cursor.fetchall(), columns=[d[0] for d in cursor.description])  "
             "# PORT: was pd.read_sql_query(sql, conn), which needs a real DB-API/SQLAlchemy connection\n   ")]),
    ],
    'src/embd/db_utils_embryos.py': [
        (DBE, 'print_rows', []),
        (DBE, 'get_columns', []),
        (DBE, 'get_all_ids', []),
        (DBE, 'create_param_names', []),
        (DBE, 'get_row', []),
        (DBE, 'insert_row', [
            ("            elif not np.isnan(colms[col]):  # if a single value (float or int)\n"
             "                names += '{pn}, '.format(pn=col)\n"
             "                values += '\\'{pv}\\', '.format(pv=colms[col])",
             "            elif not np.isnan(colms[col]):  # if a single value (float or int)\n"
             "                names += '{pn}, '.format(pn=col)\n"
             "                values += '\\'{pv}\\', '.format(pv=py2_str(colms[col]))  # PORT: py2 str() rendering")]),
        (DBE, 'update_row', [
            ("                sql += '{pn}={pv}, '.format(pn=col, pv=colms[col])",
             "                sql += '{pn}={pv}, '.format(pn=col, pv=py2_str(colms[col]))  # PORT: py2 str() rendering")]),
    ],
    'src/embd/phenoDistCalc.py': [
        (PDC, 'getPADS', []),
        (PDC, 'calc_pairwise_metrics', []),
        (PDC, 'optimizePAD', [
            ("def optimizePAD(rGroups, c):\n",
             "def optimizePAD(rGroups, c):\n"
             "    # PORT: ported as-is but already broken in the original (RNAiClass.setNDPosWeights does not exist).\n"
             "    from lmfit import Parameters, minimize, report_fit, Parameter  # PORT: was a module-level import\n")]),
        (PDC, 'get_RNAi_object_list', [
            ("def get_RNAi_object_list():", "def get_RNAi_object_list(embd_ids=range(1, 504)):"),
            ("    for i in range(1, 504):\n", "    for i in embd_ids:  # PORT: the original range(1, 504) is the default\n"),
            ("        except Exception, e:", "        except Exception as e:")]),
        (PDC, 'write_phenodistinfo_to_CSV', [
            ("def write_phenodistinfo_to_CSV(RNAi, distAll, csi):", "def write_phenodistinfo_to_CSV(RNAi, distAll, csi, fileName=None):"),
            ("    fd = open('Z:\\\\phenoDist_rank_{d}.csv'.format(d=time.strftime('%m%d%Y')), 'w')",
             "    if fileName is None:  # PORT: was 'Z:\\\\phenoDist_rank_<date>.csv'\n"
             "        fileName = FOLDER_IN + 'phenoDist_rank_{d}.csv'.format(d=time.strftime('%m%d%Y'))\n"
             "    fd = open(fileName, 'w', newline='')"),
            ("            fd.write('{0},{1},{2},{3},{4},{5},{6}\\n'.format(r1, r2, d, d20, nOv, csi[i][j], pad))",
             "            fd.write('{0},{1},{2},{3},{4},{5},{6}\\n'.format(r1, r2, py2_str(d), py2_str(d20), py2_str(nOv),\n"
             "                                                            py2_str(csi[i][j]), py2_str(pad)))  # PORT: py2 rendering (B6)"),
            ("    fd.close()", "    fd.close()\n    return fileName  # PORT")]),
        (PDC, 'write_phenodistinfo_to_mySQL', [("    from db_utils_embryos import", "    from .db_utils_embryos import")]),
        (PDC, 'run_phenodistcalc_populate', [
            ("def run_phenodistcalc_populate():", "def run_phenodistcalc_populate(fileName=None, embd_ids=range(1, 504)):"),
            ("    RNAi = get_RNAi_object_list()", "    RNAi = get_RNAi_object_list(embd_ids)"),
            ("    write_phenodistinfo_to_CSV(RNAi, distAll, csi)\n", "    fileName = write_phenodistinfo_to_CSV(RNAi, distAll, csi, fileName)\n"),
            ("    # write_phenodistinfo_to_mySQL(RNAi, distAll, csi)", "    # write_phenodistinfo_to_mySQL(RNAi, distAll, csi)\n    return fileName  # PORT")]),
    ],
    'src/embd/PADassist.py': [
        (PAD, 'getCSIFromFile', []),
        (PAD, 'getPairsParams', []),
        (PAD, 'getCommonPairs', [("        print uniqueList[i], countList[i], uniquePADvals[i]",
                                  "        print(uniqueList[i], countList[i], uniquePADvals[i])")]),
        (PAD, 'get_best_match', [("    print top_hits, seed, top_hit_vals", "    print(top_hits, seed, top_hit_vals)")]),
        (PAD, 'get_data', []),
        (PAD, 'get_data_from_mysql', [
            ("    sql = \"SELECT CONCAT('EMBD', LPAD(rnai1, 4, '0')), CONCAT('EMBD', LPAD(rnai2, 4, '0')), dist, dist2wt, dims, \" \\\n"
             "          \"CSI, PAD FROM pair_lookup\"",
             "    # PORT: the SQLite emulation has no CONCAT/LPAD; original: SELECT CONCAT('EMBD', LPAD(rnai1, 4, '0')), ...\n"
             "    sql = \"SELECT 'EMBD' || substr('0000' || rnai1, -4, 4), 'EMBD' || substr('0000' || rnai2, -4, 4), dist, \" \\\n"
             "          \"dist2wt, dims, CSI, PAD FROM pair_lookup\"")]),
        (PAD, 'makeNetwork', [("    G = nx.Graph()", "    import networkx as nx  # PORT: was a module-level import\n    G = nx.Graph()")]),
        (PAD, 'update_pad_rank', [
            ("def update_pad_rank(cursor, conn, data):", "def update_pad_rank(cursor, conn, data, embd_ids=range(1, 504)):"),
            ("    from db_utils_embryos import update_row, insert_row", "    from .db_utils_embryos import update_row, insert_row"),
            ("    for i in range(1, 504):\n", "    for i in embd_ids:  # PORT: the original range(1, 504) is the default\n")]),
        (PAD, 'output_best_matches_to_csv', PHENODIST_READS + [
            ("    savepath = 'Z://{0}_hits_matrix.csv'.format(RNAi)", "    savepath = FOLDER_IN + '{0}_hits_matrix.csv'.format(RNAi)  # PORT: was 'Z://'")]),
        (PAD, 'output_all_best_matches_to_csv_combined', [
            ("def output_all_best_matches_to_csv_combined():", "def output_all_best_matches_to_csv_combined(fileName=None, embd_ids=range(1, 504)):"),
            ("    fd = open('Z:\\\\all_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y')),\n              'w')",
             "    if fileName is None:  # PORT: was 'Z:\\\\all_best_matches_<date>.csv'\n"
             "        fileName = FOLDER_IN + 'all_best_matches_{d}.csv'.format(d=time.strftime('%m%d%Y'))\n"
             "    fd = open(fileName, 'w', newline='')"),
            ("    for i in range(1,504):\n", "    for i in embd_ids:  # PORT: the original range(1,504) is the default\n"),
            ("    fd.close()\n", "    fd.close()\n    return fileName  # PORT\n")] + PHENODIST_READS),
        (PAD, 'output_all_best_matches_to_one_csv_combined', [
            ("def output_all_best_matches_to_one_csv_combined():",
             "def output_all_best_matches_to_one_csv_combined(fileName=None, embd_ids=range(1, 504)):"),
            ("    for i in range(1,504):\n", "    for i in embd_ids:  # PORT: the original range(1,504) is the default\n"),
            ("                val = str(val)[0:5]\n", "                val = py2_str(val)[0:5]  # PORT: py2 str(float) (B6)\n", 4),
            ("    savepath = 'Z:/hits_matrix_v2_{d}.csv'.format(d=time.strftime('%m%d%Y'))\n"
             "    pd.DataFrame(all_hits_matrix).to_csv(savepath)\n",
             "    savepath = fileName if fileName is not None else \\\n"
             "        FOLDER_IN + 'hits_matrix_v2_{d}.csv'.format(d=time.strftime('%m%d%Y'))  # PORT: was 'Z:/'\n"
             "    pd.DataFrame(all_hits_matrix).to_csv(savepath)\n    return savepath  # PORT\n")] + PHENODIST_READS),
    ],
    'src/embd/RNAiClass.py': [
        (RNA, 'save_params_data', []),
        (RNA, 'load_params_data', []),
        (RNA, 'setNDPosition', [("        self.NDFlag = True\n",
                                  "        self.NDFlag = True\n        self._dist2zero = None  # PORT (perf, PERFORMANCE #1): invalidate cached getDist2Zero\n")]),
        (RNA, 'getOverlapDims', []),
        (RNA, 'getPAD', []),
        (RNA, 'get_cosine_dist', []),
        (RNA, 'get_norm_dist', []),
    ],
    'src/embd/embdFunc.py': [
        (EMB, 'get_interacting_genes', []),
        (EMB, 'check_for_interaction', [("    print [s for s in interacting if test in s]",
                                         "    print([s for s in interacting if test in s])")]),
        (EMB, 'getEMBDNumber', [("    elif ind != [None]:",
                                 "    elif ind.size == 1:  # PORT: py2 `ind != [None]` is truthy for a 1-element array, falsy if empty"),
                                ("    data = readCSV(FILENAME_GENE_NAMES, ',')[1:]\n    dataA = np.array(data)\n",
                                 "    data, dataA = _gene_key_rows()  # PORT (perf): was readCSV(FILENAME_GENE_NAMES, ',')[1:] and np.array(data)\n")]),
    ],
    'src/embd/myMath.py': [
        (MTH, 'nan_helper', []),
    ],
}


def find_block(src, name):
    """(start, end) character offsets of the last ``def name(`` block (with its indentation) in src."""
    lines = src.splitlines(keepends=True)
    starts = [i for i, l in enumerate(lines) if re.match(r'^\s*def %s\(' % re.escape(name), l)]
    if not starts:
        raise KeyError(name)
    s = starts[-1]
    ind = len(lines[s]) - len(lines[s].lstrip())
    e = s + 1
    while e < len(lines):
        l = lines[e]
        if l.strip() and (len(l) - len(l.lstrip())) <= ind and not l.lstrip().startswith('#'):
            break
        e += 1
    while e > s + 1 and not lines[e - 1].strip():  # keep trailing blank lines outside the block
        e -= 1
    off = [0]
    for l in lines:
        off.append(off[-1] + len(l))
    return off[s], off[e]


def build(orig_path, name, edits):
    src = open(os.path.join(REPO, orig_path), encoding='latin-1').read()
    a, b = find_block(src, name)
    body = src[a:b]
    for ed in edits:
        old, new = ed[0], ed[1]
        n = ed[2] if len(ed) > 2 else 1
        c = body.count(old)
        if c != n:
            raise AssertionError('%s:%s edit anchor found %d times (expected %d): %r' % (orig_path, name, c, n, old[:70]))
        body = body.replace(old, new)
    if not body.endswith('\n'):
        body += '\n'
    return body


def main():
    for port_path, entries in SPEC.items():
        path = os.path.join(REPO, port_path)
        src = open(path, encoding='utf-8').read()
        for orig_path, name, edits in entries:
            body = build(orig_path, name, edits)
            a, b = find_block(src, name)
            src = src[:a] + body + src[b:]
        with open(path, 'w', encoding='utf-8') as f:
            f.write(src)
        print('rebuilt %d functions in %s' % (len(entries), port_path))


if __name__ == '__main__':
    main()

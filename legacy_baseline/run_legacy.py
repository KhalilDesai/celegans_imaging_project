# -*- coding: utf-8 -*-
"""Driver that runs the ORIGINAL, unmodified Python 2.7 downstream code (legacy baseline).

Runs inside docker/legacy.Dockerfile. It never edits the original sources. It only
  * puts import stubs (legacy_baseline/shims) and the data layer (src/) ahead of them on sys.path,
  * supplies names the original code expects from the lab environment but that are absent from
    the public repo (``varLookup.MYSQL_DB_CREDENTIALS``, ``db_utils_embryos.mySQLpath/login``,
    ``fileLookup.SKIPPED``, the module-global ``db_utils_genes.conn``),
  * selects the parameter set the way the authors did, by swapping ``params_use_assign.get_params_use``
    before ``RNAiClass`` is imported (one process per parameter set),
  * renames date-stamped outputs to the file names the downstream readers expect,
  * runs the original functions in the order reconstructed in docs/PORTING_PLAN.md section 2.

Every deviation from "call the original function as-is" is listed in docs/DECISIONS_AND_ISSUES.md
(section L) and marked ``LEGACY-DRIVER`` below.

Usage (inside the container, cwd = work dir that contains the ``Z:`` directory):
    python run_legacy.py init-db
    python run_legacy.py phenodist {opt,gls,ms,all}
    python run_legacy.py vectors {opt,gls,ms,all}
    python run_legacy.py bestmatches
    python run_legacy.py padrank
    python run_legacy.py derived
    python run_legacy.py matrices {opt,gls,ms,all}
"""
from __future__ import print_function

import csv
import glob
import json
import os
import shutil
import sqlite3
import sys
import time

REPO = os.environ.get('EMBD_REPO', '/repo')
WORK = os.environ.get('EMBD_WORK', '/work')
Z = 'Z:'  # FOLDER_IN = 'Z:/' in varLookup resolves to <cwd>/Z:/ on Linux

PARAMSET_FUNCS = {'opt': '_group_rank_opt', 'gls': '_gls_params', 'ms': '_ms_params', 'all': '_all'}
# File names the original readers expect (fileLookup.FILENAME_PAD, PADassist.output_all_best_matches_*).
PHENODIST_NAMES = {'opt': 'phenoDist_rank_04142022.csv',
                   'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
                   'ms': 'phenoDist_rank_MS_ONLY05102022.csv',
                   'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}
UNIVERSE = list(range(1, 504))  # DECISIONS A10


def log_timing(step, seconds, extra=None):
    path = os.path.join(Z, 'timings_legacy.jsonl')
    rec = {'step': step, 'seconds': round(seconds, 3)}
    rec.update(extra or {})
    with open(path, 'a') as f:
        f.write(json.dumps(rec, sort_keys=True) + '\n')
    print('[timing]', json.dumps(rec, sort_keys=True))


def bootstrap(paramset=None):
    shims = os.path.join(REPO, 'legacy_baseline', 'shims')
    for p in reversed([shims, os.path.join(REPO, 'src'),
                       os.path.join(REPO, 'development'), os.path.join(REPO, 'embryos_analyze')]):
        if p not in sys.path:
            sys.path.insert(0, p)
    os.chdir(WORK)
    for d in (Z, os.path.join(Z, 'LOG')):
        if not os.path.isdir(d):
            os.makedirs(d)
    import varLookup
    varLookup.MYSQL_DB_CREDENTIALS = None  # LEGACY-DRIVER: commented out in the public varLookup.py
    import db_utils_embryos
    db_utils_embryos.mySQLpath = 'sqlite-emulation'  # LEGACY-DRIVER: undefined in the public repo
    db_utils_embryos.login = 'embd'
    import fileLookup
    if not hasattr(fileLookup, 'SKIPPED'):
        fileLookup.SKIPPED = []  # LEGACY-DRIVER: undefined in the public repo (DECISIONS C5/G6, Q4)
    import params_use_assign
    if paramset is not None:
        # LEGACY-DRIVER: the authors switched parameter sets by editing get_params_use(); same effect.
        params_use_assign.get_params_use = getattr(params_use_assign, PARAMSET_FUNCS[paramset])
    if 'RNAiClass' in sys.modules:
        raise RuntimeError('RNAiClass imported before parameter set was selected')


def _move_single(pattern, dest):
    found = sorted(f for f in glob.glob(pattern) if os.path.abspath(f) != os.path.abspath(dest))
    if len(found) != 1:
        raise RuntimeError('expected exactly one file matching %r, found %r' % (pattern, found))
    shutil.move(found[0], dest)
    return dest


def cmd_init_db():
    sys.path.insert(0, os.path.join(REPO, 'src'))
    from embd_data.mysql_emulation import init_db
    os.chdir(WORK)
    if not os.path.isdir(Z):
        os.makedirs(Z)
    types = init_db(os.environ['EMBD_SQLITE'], os.environ['EMBD_GENES_CSV'])
    shutil.copy(os.environ['EMBD_GENE_KEY_CSV'], os.path.join(Z, 'gene_key.csv'))
    print('genes table: %d columns (%d FLOAT)' % (len(types), sum(1 for t in types.values() if t == 'FLOAT')))


def cmd_phenodist(paramset):
    bootstrap(paramset)
    import phenoDistCalc
    t0 = time.time()
    phenoDistCalc.run_phenodistcalc_populate()
    dt = time.time() - t0
    # write_phenodistinfo_to_CSV opens 'Z:\\phenoDist_rank_<date>.csv' -> a file literally named so in cwd
    dest = _move_single('Z:\\phenoDist_rank_*.csv', os.path.join(Z, PHENODIST_NAMES[paramset]))
    log_timing('phenodist_' + paramset, dt, {'output': dest})


def cmd_vectors(paramset):
    """Instrumentation (not part of the original pipeline): dump per-gene ND positions."""
    bootstrap(paramset)
    import numpy as np
    from RNAiClass import RNAiClass, PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS
    t0 = time.time()
    pos, d20 = [], []
    for i in UNIVERSE:
        r = RNAiClass(i, verb=False)
        r.setNDPosition()
        pos.append(r.posND)
        d20.append(r.getDist2Zero())
    np.save(os.path.join(Z, 'posND_%s.npy' % paramset), np.array(pos))
    np.save(os.path.join(Z, 'dist2zero_%s.npy' % paramset), np.array(d20))
    with open(os.path.join(Z, 'params_used_%s.json' % paramset), 'w') as f:
        json.dump({'GLS': list(PARAM_NAMES_USE_GLS), 'MS': list(PARAM_NAMES_USE_MS)}, f, indent=1)
    log_timing('vectors_' + paramset, time.time() - t0)


def cmd_bestmatches():
    bootstrap('opt')
    import PADassist
    if os.environ.get('EMBD_LEGACY_CACHE_READCSV', '1') == '1':
        # LEGACY-DRIVER (performance only): the original re-reads four ~253k-row CSVs inside a
        # 503-iteration loop. readCSV is a pure function of an unchanging file, and its result is
        # only read (PADFilter copies it with np.array), so memoising it cannot change outputs.
        _orig, _cache = PADassist.readCSV, {}

        def cached_readCSV(fileName, delimiter=' '):
            key = (fileName, delimiter)
            if key not in _cache:
                _cache[key] = _orig(fileName, delimiter)
            return _cache[key]
        PADassist.readCSV = cached_readCSV
    t0 = time.time()
    PADassist.output_all_best_matches_to_csv_combined()
    log_timing('bestmatches_combined', time.time() - t0)
    _move_single('Z:\\all_best_matches_*.csv', os.path.join(Z, 'all_best_matches_legacy.csv'))
    t0 = time.time()
    PADassist.output_all_best_matches_to_one_csv_combined()
    log_timing('bestmatches_one_csv', time.time() - t0)
    _move_single(os.path.join(Z, 'hits_matrix_v2_*.csv'), os.path.join(Z, 'hits_matrix_v2_legacy.csv'))


def _export_table(table, dest, columns=None):
    db = sqlite3.connect(os.environ['EMBD_SQLITE'])
    cur = db.execute('SELECT %s FROM %s' % (', '.join('`%s`' % c for c in columns) if columns else '*', table))
    names = [d[0] for d in cur.description]
    with open(dest, 'wb') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(names)
        for r in cur:
            w.writerow([(repr(v) if isinstance(v, float) else (v.encode('utf-8') if isinstance(v, unicode) else v))  # noqa: F821
                        for v in r])
    db.close()


def cmd_padrank():
    bootstrap('opt')
    import PADassist
    from db_utils_embryos import initialize_db
    conn, curs = initialize_db()
    data = PADassist.get_data()  # get_data_from_csv(): FILENAME_PAD = Z:/phenoDist_rank_04142022.csv
    t0 = time.time()
    PADassist.update_pad_rank(curs, conn, data)
    log_timing('padrank', time.time() - t0)
    _export_table('pad_rank', os.path.join(Z, 'pad_rank_legacy.csv'))


def _update_top_hits_PAD_all(cursor, conn, embd_ids):
    """LEGACY-DRIVER: verbatim body of db_utils_genes.update_top_hits_PAD; only the hand-edited loop
    list ([1406, 674, ...] in the public file) is replaced by ``embd_ids``."""
    import numpy as np
    from embdFunc import getGeneCommonName
    thr = 0.5
    for i in embd_ids:
        print('starting {0}'.format(i))
        seed = i
        ranks_considered = 'rank1, rank2, rank3, rank4, rank5, rank6, rank7, rank8, rank9, rank10'
        pad_considered = 'pad1, pad2, pad3, pad4, pad5, pad6, pad7, pad8, pad9, pad10'
        sql = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(ranks_considered, seed)
        sql2 = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(pad_considered, seed)
        cursor.execute(sql)
        hits = cursor.fetchall()
        hits = np.array(hits[0])
        cursor.execute(sql2)
        pad = cursor.fetchall()
        pad = np.array(pad[0])
        high_pad = pad > thr
        inds = np.where(high_pad == True)  # noqa: E712
        ranks_use = hits[inds]
        genes = []
        for item in ranks_use:
            genes.append(getGeneCommonName(item))
        genes_use = ', '.join(genes)
        sql = "SELECT id, embd_id, top_hits_PAD FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, top_hits_PAD = cursor.fetchall()[0]
        sql2 = "UPDATE genes SET top_hits_PAD=\"{h}\" WHERE id={id};".format(h=genes_use, id=id)
        cursor.execute(sql2)
        conn.commit()


def _update_column_reciprocal_matches_all(cursor, conn, embd_ids):
    """LEGACY-DRIVER: verbatim body of db_utils_genes.update_column_reciprocal_matches; only the
    hand-edited ``range(310, 311)`` is replaced by ``embd_ids``."""
    from embdFunc import getGeneCommonName
    from db_utils_genes import get_reciprocal_matches
    for i in embd_ids:
        best_hits = get_reciprocal_matches(i, 10, cursor, thr=0.55)
        top_five = best_hits[0:10]
        genes = []
        for item in best_hits:
            genes.append(getGeneCommonName(item))
        top_five_genes = ', '.join(genes[0:10])
        print((i, top_five, top_five_genes))
        sql = "SELECT id, embd_id, movement_GLS, movement_MS, WT_like_or_phenotype, reciprocal_matches FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, movementG, movementM, WT_like_or_phenotype, reciprocal_matches = cursor.fetchall()[0]
        if movementG == 0 or movementM == 0:
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1 and "Dev. Delay" in WT_like_or_phenotype:
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        else:
            wt_note = 'Embryos do not exhibit an arrest phenotype'
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=wt_note, id=id)
            cursor.execute(sql2)
            conn.commit()


def cmd_derived():
    bootstrap('opt')
    import db_utils_genes
    from db_utils_embryos import initialize_db
    conn, curs = initialize_db()
    db_utils_genes.conn = conn  # LEGACY-DRIVER: the update_* functions use the __main__-level global
    # Order (DECISIONS Q5, PORTING_PLAN section 2): auto_description -> WT_like -> top hits -> reciprocal.
    t0 = time.time()
    db_utils_genes.update_automated_description()  # -> groupAssesment.get_all_closest_groups, EMBD 37..503
    log_timing('derived_auto_description', time.time() - t0)
    t0 = time.time()
    db_utils_genes.update_WT_like_or_phenotype(curs)  # original loop is already range(1, 504)
    _update_top_hits_PAD_all(curs, conn, UNIVERSE)
    _update_column_reciprocal_matches_all(curs, conn, UNIVERSE)
    log_timing('derived_hits_recip', time.time() - t0)
    conn.close()
    _export_table('genes', os.path.join(Z, 'genes_derived_legacy.csv'),
                  ['embd_id', 'movement_GLS', 'movement_MS', 'auto_description', 'WT_like_or_phenotype',
                   'top_hits_PAD', 'reciprocal_matches'])


def cmd_matrices(paramset):
    bootstrap(paramset)
    import groupAssesment
    t0 = time.time()
    groupAssesment.make_UMAP_content_files()
    log_timing('matrices_' + paramset, time.time() - t0)
    for kind in ('PAD_Matrix', 'Dist_Matrix', 'Matrix_labels', 'group_Lists', 'group_Labels', 'group_colors',
                 'group_Labels_toPlot'):
        shutil.move(os.path.join(Z, '%s_msParams.npy' % kind), os.path.join(Z, '%s_%s.npy' % (kind, paramset)))


def main(argv):
    cmd = argv[1]
    if cmd == 'init-db':
        cmd_init_db()
    elif cmd in ('phenodist', 'vectors', 'matrices'):
        if len(argv) < 3 or argv[2] not in PARAMSET_FUNCS:
            raise SystemExit('usage: %s %s {opt,gls,ms,all}' % (argv[0], cmd))
        {'phenodist': cmd_phenodist, 'vectors': cmd_vectors, 'matrices': cmd_matrices}[cmd](argv[2])
    elif cmd == 'bestmatches':
        cmd_bestmatches()
    elif cmd == 'padrank':
        cmd_padrank()
    elif cmd == 'derived':
        cmd_derived()
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    main(sys.argv)

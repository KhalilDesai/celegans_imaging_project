"""Run the ported downstream pipeline (Python 3), mirroring legacy_baseline/run_legacy.py step by step.

    python -m embd.pipeline init-db  --genes-csv build/canonical/genes_table.csv --gene-key build/canonical/gene_key.csv
    python -m embd.pipeline phenodist {opt,gls,ms,all}
    python -m embd.pipeline vectors   {opt,gls,ms,all}
    python -m embd.pipeline matrices  {opt,gls,ms,all}
    python -m embd.pipeline bestmatches
    python -m embd.pipeline padrank
    python -m embd.pipeline derived
    python -m embd.pipeline umap      {opt,gls,ms,all}
    python -m embd.pipeline all       # every step, parameter sets in separate processes

Environment: EMBD_OUT_DIR (outputs; the original 'Z:/'), EMBD_SQLITE (emulated MySQL file).
Each parameter-set step runs in its own process with EMBD_PARAMSET set, because the original
modules fix the parameter set at import time (DECISIONS L5).
Output file names are those of the legacy baseline, so outputs can be diffed with tests/golden.
"""
import json
import os
import shutil
import subprocess
import sys
import time

SETS = ('opt', 'gls', 'ms', 'all')
PHENODIST = {'opt': 'phenoDist_rank_04142022.csv', 'gls': 'phenoDist_rank_GLS_ONLY05102022.csv',
             'ms': 'phenoDist_rank_MS_ONLY05102022.csv', 'all': 'phenoDist_rank_ALL_PARAMS05102022.csv'}
UNIVERSE = range(1, 504)  # DECISIONS A10


def _out():
    from embd import config
    return config.folder_in()


def _timing(step, seconds):
    with open(os.path.join(_out(), 'timings_port.jsonl'), 'a') as f:
        f.write(json.dumps({'step': step, 'seconds': round(seconds, 3)}) + '\n')


def _require_paramset(ps):
    if os.environ.get('EMBD_PARAMSET', 'opt') != ps:
        raise RuntimeError('run this step in a process with EMBD_PARAMSET=%s (use `pipeline all` or set it)' % ps)


def init_db(genes_csv, gene_key):
    from embd_data.mysql_emulation import init_db as _init
    _init(os.environ['EMBD_SQLITE'], genes_csv)
    shutil.copy(gene_key, os.path.join(_out(), 'gene_key.csv'))


def phenodist(ps, embd_ids=UNIVERSE):
    _require_paramset(ps)
    from embd import phenoDistCalc
    t0 = time.time()
    phenoDistCalc.run_phenodistcalc_populate(fileName=_out() + PHENODIST[ps], embd_ids=embd_ids)
    _timing('phenodist_' + ps, time.time() - t0)


def vectors(ps, embd_ids=UNIVERSE):
    _require_paramset(ps)
    import numpy as np
    from embd.RNAiClass import PARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS, RNAiClass
    pos, d20 = [], []
    for i in embd_ids:
        r = RNAiClass(i, verb=False)
        r.setNDPosition()
        pos.append(r.posND)
        d20.append(r.getDist2Zero())
    np.save(_out() + 'posND_%s.npy' % ps, np.array(pos))
    np.save(_out() + 'dist2zero_%s.npy' % ps, np.array(d20))
    with open(_out() + 'params_used_%s.json' % ps, 'w') as f:
        json.dump({'GLS': list(PARAM_NAMES_USE_GLS), 'MS': list(PARAM_NAMES_USE_MS)}, f, indent=1)


def matrices(ps):
    _require_paramset(ps)
    from embd import groupAssesment
    t0 = time.time()
    groupAssesment.make_UMAP_content_files()
    _timing('matrices_' + ps, time.time() - t0)
    for kind in ('PAD_Matrix', 'Dist_Matrix', 'Matrix_labels', 'group_Lists', 'group_Labels', 'group_colors',
                 'group_Labels_toPlot'):
        os.replace(_out() + '%s_msParams.npy' % kind, _out() + '%s_%s.npy' % (kind, ps))


def bestmatches(embd_ids=UNIVERSE):
    from embd import PADassist
    t0 = time.time()
    PADassist.output_all_best_matches_to_csv_combined(_out() + 'all_best_matches_port.csv', embd_ids)
    _timing('bestmatches_combined', time.time() - t0)
    t0 = time.time()
    PADassist.output_all_best_matches_to_one_csv_combined(_out() + 'hits_matrix_v2_port.csv', embd_ids)
    _timing('bestmatches_one_csv', time.time() - t0)


def padrank(embd_ids=UNIVERSE):
    from embd import PADassist
    from embd.db_utils_embryos import initialize_db
    conn, curs = initialize_db()
    t0 = time.time()
    PADassist.update_pad_rank(curs, conn, PADassist.get_data(), embd_ids)
    _timing('padrank', time.time() - t0)
    export_table('pad_rank', _out() + 'pad_rank_port.csv')


def derived(embd_ids=UNIVERSE):
    from embd import db_utils_genes
    from embd.db_utils_embryos import initialize_db
    conn, curs = initialize_db()
    t0 = time.time()
    db_utils_genes.update_automated_description()  # groupAssesment.get_all_closest_groups, EMBD 37..503
    db_utils_genes.update_WT_like_or_phenotype(curs, conn)
    db_utils_genes.update_top_hits_PAD(curs, conn, embd_ids=embd_ids)
    db_utils_genes.update_column_reciprocal_matches(curs, conn, embd_ids=embd_ids)
    _timing('derived', time.time() - t0)
    conn.close()
    export_table('genes', _out() + 'genes_derived_port.csv',
                 ['embd_id', 'movement_GLS', 'movement_MS', 'auto_description', 'WT_like_or_phenotype',
                  'top_hits_PAD', 'reciprocal_matches'])


def umap_step(ps):
    import numpy as np
    from embd import makeUMAPs_V2
    from embd.embdFunc import getGeneCommonName
    data = np.load(_out() + 'PAD_Matrix_%s.npy' % ps)
    t0 = time.time()
    pts = makeUMAPs_V2.run(data, query=1, common_name=getGeneCommonName(1),
                           svg_path=_out() + 'UMAP_PAD_%s_EMBD0001.svg' % ps)
    _timing('umap_' + ps, time.time() - t0)
    np.save(_out() + 'umap_points_%s.npy' % ps, pts)


def export_table(table, dest, columns=None):
    import csv
    import sqlite3
    db = sqlite3.connect(os.environ['EMBD_SQLITE'])
    cur = db.execute('SELECT %s FROM %s' % (', '.join('`%s`' % c for c in columns) if columns else '*', table))
    names = [d[0] for d in cur.description]
    with open(dest, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(names)
        for r in cur:
            w.writerow([repr(v) if isinstance(v, float) else v for v in r])
    db.close()


def _sub(*args, paramset=None):
    env = dict(os.environ)
    # Do not rely on the editable-install .pth file: macOS can flag it hidden and Python >= 3.12.8
    # then skips it (DECISIONS P2-8).
    src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env['PYTHONPATH'] = src + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    if paramset:
        env['EMBD_PARAMSET'] = paramset
    subprocess.check_call([sys.executable, '-m', 'embd.pipeline'] + list(args), env=env)


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    cmd, rest = argv[1], argv[2:]
    if cmd == 'init-db':
        init_db(rest[rest.index('--genes-csv') + 1], rest[rest.index('--gene-key') + 1])
    elif cmd in ('phenodist', 'vectors', 'matrices'):
        {'phenodist': phenodist, 'vectors': vectors, 'matrices': matrices}[cmd](rest[0])
    elif cmd == 'umap':
        umap_step(rest[0])
    elif cmd == 'bestmatches':
        bestmatches()
    elif cmd == 'padrank':
        padrank()
    elif cmd == 'derived':
        derived()
    elif cmd == 'all':
        for s in SETS:
            _sub('phenodist', s, paramset=s)
        _sub('bestmatches')
        _sub('padrank')
        _sub('derived')
        for s in SETS:
            _sub('vectors', s, paramset=s)
            _sub('matrices', s, paramset=s)
            _sub('umap', s)
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    main(sys.argv)

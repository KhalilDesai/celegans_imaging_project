"""Tests for the spreadsheet -> canonical CSV -> MySQL-emulation data layer."""
import csv
import math
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from embd_data import mysql_emulation as me  # noqa: E402
from embd_data import sheet  # noqa: E402

SHEET = os.path.join(REPO, 'data', 'All_genes_mySQL_table.xlsx')
INDEX = os.path.join(REPO, 'data', 'Index_EMBD_ID_Gene_names.xlsx')


@pytest.fixture(scope='module')
def canonical(tmp_path_factory):
    pytest.importorskip('openpyxl')
    out = tmp_path_factory.mktemp('canon')
    sheet.convert_genes_table(SHEET, str(out / 'genes_table.csv'))
    sheet.write_gene_key(INDEX, str(out / 'gene_key.csv'))
    return out


@pytest.fixture()
def db(canonical, tmp_path, monkeypatch):
    path = str(tmp_path / 'db.sqlite')
    me.init_db(path, str(canonical / 'genes_table.csv'))
    monkeypatch.setenv('EMBD_SQLITE', path)
    return path


def test_canonical_csv_round_trips_every_float_exactly(canonical):
    header, body = sheet._read_xlsx(SHEET)
    with open(canonical / 'genes_table.csv', newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f))
    out_h = rows[0]
    by_embd = {int(r[out_h.index('embd_id')]): r for r in rows[1:]}
    checked = 0
    for r in body:
        rec = dict(zip(header, r))
        out = by_embd[int(rec['embd_id'])]
        for name, v in rec.items():
            if name is None:
                continue
            got = out[out_h.index(name)]
            if v is None or v == 'NULL':
                assert got == ''
            elif isinstance(v, float):
                assert float(got) == v  # bit-exact
                checked += 1
            elif isinstance(v, int):
                assert int(got) == v
    assert checked == 60764  # every float cell in the checksum-pinned sheet


def test_control_row_gets_hardcoded_mysql_id(canonical):
    with open(canonical / 'genes_table.csv', newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    ctrl = [r for r in rows if r['embd_id'] == '0']
    assert len(ctrl) == 1 and ctrl[0]['id'] == '2627' and ctrl[0]['version'] == sheet.DATA_VERSION
    assert len(rows) == 514


def test_translate_mysql_literals():
    assert me.translate_mysql('UPDATE genes SET a="x\'y" WHERE id=1') == "UPDATE genes SET a='x''y' WHERE id=1"
    assert me.translate_mysql("SELECT 'a\\nb'") == "SELECT 'a\nb'"
    assert me.translate_mysql('SET a=""') == "SET a=''"


def test_mysql_float_render():
    assert me.mysql_float_render(0.856809720320) == 0.85681
    assert me.mysql_float_render(125251000.4) == 125251000.0
    assert me.mysql_float_render(24.7983) == 24.7983


def test_get_row_like_original(db):
    conn = me.connect()
    cur = conn.cursor()
    cur.execute('SHOW COLUMNS FROM genes')
    cols = [c[0] for c in cur.fetchall()]
    assert 'aG_GLS' in cols and 'version' in cols and 'bY_GLS' not in cols
    cur.execute('SELECT id FROM genes WHERE embd_id=1')
    (gid,), = cur.fetchall()
    cur.execute('SELECT aG_GLS, aG_MS, version FROM genes WHERE id={}'.format(gid))
    row = cur.fetchall()[0]
    assert [d[0] for d in cur.description] == ['aG_GLS', 'aG_MS', 'version']
    assert row[2] == sheet.DATA_VERSION


def test_pad_rank_float_emulation(db):
    conn = me.connect()
    cur = conn.cursor()
    cur.execute("INSERT INTO pad_rank (rnai, rank1, pad1) VALUES ('5', '76', '0.856809720320')")
    conn.commit()
    cur.execute('SELECT rank1, pad1 FROM PAD_RANK WHERE rnai=5;')
    assert cur.fetchall() == [(76, 0.85681)]


def test_unsupported_sql_fails_loudly(db):
    cur = me.connect().cursor()
    with pytest.raises(me.Error):
        cur.execute('SELECT CONCAT_WS(1) FROM nonexistent_table')


def test_gene_key_layout(canonical):
    with open(canonical / 'gene_key.csv', newline='', encoding='utf-8') as f:
        rows = list(csv.reader(f))
    assert rows[0] == ['EMBD_ID', 'Sequence_Name', 'Common_Name']
    assert rows[1] == ['EMB_P0001', 'K02F2.3', 'teg-4']
    assert len(rows) == 1 + 513
    assert not any(math.isnan(float(r[0][5:])) for r in rows[1:])

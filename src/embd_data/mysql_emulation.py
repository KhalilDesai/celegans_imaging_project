# -*- coding: utf-8 -*-
"""A ``MySQLdb``-compatible stand-in backed by a local SQLite file.

Lets the *unmodified* original code (``db_utils_embryos``, ``db_utils_genes``, ``PADassist``,
``phenoDistCalc`` ...) run without a MySQL server. Compatible with Python 2.7 and Python 3.

Faithfulness notes (see docs/DECISIONS_AND_ISSUES.md, section "D"):

* Only the SQL shapes used by the in-scope code are supported. Anything SQLite cannot parse
  raises immediately, so an unsupported path cannot fail silently.
* MySQL string literals ('...' and "..." with backslash escapes, the default ``sql_mode``) are
  rewritten as SQLite single-quoted literals.
* ``SHOW COLUMNS FROM t`` is emulated with ``PRAGMA table_info``.
* MySQL ``FLOAT`` columns are single precision, and this server rendered them with 6
  significant digits (that rendering is visible in the published dump, DECISIONS A6). Reading
  a FLOAT column can therefore be emulated as ``float('%.6g' % float32(x))``. The mode is set
  per table via ``EMBD_MYSQL_FLOAT_MODE_<TABLE>`` = ``asis`` | ``mysql``. Defaults:
  ``genes=asis`` (the sheet already holds the rendered values), ``pad_rank=mysql``,
  ``pair_lookup=mysql`` (values written by the pipeline itself).
* Writes go to the SQLite file named by ``EMBD_SQLITE`` and persist across connections, as a
  database server would. The spreadsheet itself is never modified.
"""
from __future__ import print_function

import csv
import io
import os
import re
import sqlite3
import struct
import sys

PY2 = sys.version_info[0] == 2
if PY2:  # pragma: no cover - exercised only in the legacy container
    text_type = unicode  # noqa: F821
else:
    text_type = str

DEFAULT_FLOAT_MODES = {'genes': 'asis', 'pad_rank': 'mysql', 'pair_lookup': 'mysql'}

# MySQL schema of the derived tables, reconstructed from db_utils_genes.create_columns_pad_rank
# and phenoDistCalc.write_phenodistinfo_to_mySQL.
PAD_RANK_DDL = ('CREATE TABLE pad_rank (id INTEGER PRIMARY KEY AUTOINCREMENT, rnai INTEGER, ' +
                ', '.join('rank{0} INTEGER, pad{0} FLOAT'.format(i) for i in range(1, 51)) + ')')
PAIR_LOOKUP_DDL = ('CREATE TABLE pair_lookup (id INTEGER PRIMARY KEY AUTOINCREMENT, rnai1 INTEGER, '
                   'rnai2 INTEGER, dist FLOAT, dist2wt FLOAT, dims INTEGER, CSI FLOAT, PAD FLOAT)')


class Error(Exception):
    pass


OperationalError = ProgrammingError = Error


def mysql_float_render(x):
    """Value MySQLdb returns for a MySQL FLOAT column holding ``x`` (float32, '%.6g' text)."""
    if x is None:
        return None
    f32 = struct.unpack('f', struct.pack('f', float(x)))[0]
    return float('%.6g' % f32)


def _float_mode(table):
    t = table.lower()
    return os.environ.get('EMBD_MYSQL_FLOAT_MODE_' + t.upper(), DEFAULT_FLOAT_MODES.get(t, 'asis'))


_SHOW_COLUMNS = re.compile(r'^\s*SHOW\s+COLUMNS\s+FROM\s+`?(\w+)`?\s*;?\s*$', re.I)
_FROM_TABLE = re.compile(r'\bFROM\s+`?(\w+)`?', re.I)


def translate_mysql(sql):
    """Rewrite MySQL string literals ('...' or "..." with backslash escapes) as SQLite literals."""
    out, i, n = [], 0, len(sql)
    escapes = {'0': '\0', 'b': '\b', 'n': '\n', 'r': '\r', 't': '\t', 'Z': '\x1a'}
    while i < n:
        c = sql[i]
        if c in ('"', "'"):
            q, i, buf = c, i + 1, []
            while True:
                if i >= n:
                    raise Error('unterminated string literal in: %r' % sql)
                ch = sql[i]
                if ch == '\\' and i + 1 < n:
                    buf.append(escapes.get(sql[i + 1], sql[i + 1]))
                    i += 2
                    continue
                if ch == q:
                    if i + 1 < n and sql[i + 1] == q:  # doubled quote
                        buf.append(q)
                        i += 2
                        continue
                    i += 1
                    break
                buf.append(ch)
                i += 1
            out.append("'" + ''.join(buf).replace("'", "''") + "'")
        else:
            out.append(c)
            i += 1
    return ''.join(out)


class Cursor(object):
    def __init__(self, conn):
        self._conn = conn
        self._cur = conn._db.cursor()
        self._rows = None
        self.description = None
        self.rowcount = -1

    def execute(self, sql, args=None):
        if args:
            raise Error('parameterised queries are not used by the original code: %r' % sql)
        if PY2 and isinstance(sql, str):
            sql = sql.decode('utf-8')
        m = _SHOW_COLUMNS.match(sql)
        if m:
            info = self._conn._table_info(m.group(1))
            if not info:
                raise Error("Table '%s' doesn't exist" % m.group(1))
            self._rows = [(self._s(name), self._s(decl), 'YES', '', None, '') for name, decl in info]
            self.description = tuple((d, None, None, None, None, None, None)
                                     for d in ('Field', 'Type', 'Null', 'Key', 'Default', 'Extra'))
            self.rowcount = len(self._rows)
            return self.rowcount
        try:
            self._cur.execute(translate_mysql(sql))
        except sqlite3.Error as e:
            raise Error('%s (SQL: %s)' % (e, sql))
        self.description = self._cur.description
        self.rowcount = self._cur.rowcount
        if self.description is None:
            self._rows = None
            return self.rowcount
        rows = self._cur.fetchall()
        mt = _FROM_TABLE.search(sql)
        if mt and _float_mode(mt.group(1)) == 'mysql':
            types = dict((name.lower(), decl.upper()) for name, decl in self._conn._table_info(mt.group(1)))
            fcols = [j for j, d in enumerate(self.description) if types.get(d[0].lower()) == 'FLOAT']
            if fcols:
                rows = [tuple(mysql_float_render(v) if j in fcols and v is not None else v
                              for j, v in enumerate(r)) for r in rows]
        self._rows = [tuple(self._s(v) for v in r) for r in rows]
        self.rowcount = len(self._rows)
        return self.rowcount

    @staticmethod
    def _s(v):
        if PY2 and isinstance(v, text_type):
            return v.encode('utf-8')  # MySQLdb on Python 2 returns byte strings
        return v

    def fetchall(self):
        rows, self._rows = self._rows or [], []
        return rows

    def fetchone(self):
        if not self._rows:
            return None
        return self._rows.pop(0)

    def close(self):
        self._cur.close()

    def __iter__(self):
        return iter(self.fetchall())


class Connection(object):
    def __init__(self, path):
        self._db = sqlite3.connect(path)
        self._info_cache = {}

    def _table_info(self, table):
        key = table.lower()
        if key not in self._info_cache:
            rows = self._db.execute('PRAGMA table_info(`%s`)' % table).fetchall()
            self._info_cache[key] = [(r[1], r[2]) for r in rows]
        return self._info_cache[key]

    def cursor(self):
        return Cursor(self)

    def commit(self):
        self._db.commit()

    def rollback(self):
        self._db.rollback()

    def close(self):
        self._db.close()


def connect(*args, **kwargs):
    """Drop-in for ``MySQLdb.connect``. Connection arguments are ignored; the database is ``$EMBD_SQLITE``."""
    path = os.environ.get('EMBD_SQLITE')
    if not path or not os.path.exists(path):
        raise Error('EMBD_SQLITE must point to a database created by embd_data.mysql_emulation.init_db')
    return Connection(path)


def _read_csv_rows(path):
    if PY2:
        with open(path, 'rb') as f:
            rows = list(csv.reader(f))
        return [[c.decode('utf-8') for c in r] for r in rows]
    with io.open(path, 'r', encoding='utf-8', newline='') as f:
        return list(csv.reader(f))


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def init_db(sqlite_path, genes_csv):
    """Create the SQLite database from the canonical genes CSV (see embd_data.sheet)."""
    if os.path.exists(sqlite_path):
        os.remove(sqlite_path)
    rows = _read_csv_rows(genes_csv)
    header, body = rows[0], rows[1:]
    decl = []
    for j, name in enumerate(header):
        col = [r[j] for r in body if r[j] != '']
        if name in ('id', 'embd_id'):
            decl.append('INTEGER')
        elif name == 'version':
            decl.append('TEXT')
        elif all(_is_number(v) for v in col):
            decl.append('FLOAT')
        else:
            decl.append('TEXT')
    db = sqlite3.connect(sqlite_path)
    cols_sql = ', '.join('`%s` %s%s' % (n, t, ' PRIMARY KEY' if n == 'id' else '') for n, t in zip(header, decl))
    db.execute('CREATE TABLE genes (%s)' % cols_sql)

    def conv(v, t):
        if v == '':
            return None
        if t == 'INTEGER':
            return int(v)
        if t == 'FLOAT':
            return float(v)
        return v
    db.executemany('INSERT INTO genes VALUES (%s)' % ', '.join('?' * len(header)),
                   [tuple(conv(v, t) for v, t in zip(r, decl)) for r in body])
    db.execute(PAD_RANK_DDL)
    db.execute(PAIR_LOOKUP_DDL)
    db.commit()
    db.close()
    return dict(zip(header, decl))


if __name__ == '__main__':
    init_db(sys.argv[1], sys.argv[2])
    print('created', sys.argv[1])

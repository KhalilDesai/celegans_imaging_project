"""
Port of embryos_analyze/db_utils_embryos.py: generic row access used by the downstream code.

``initialize_db`` connects to the SQLite-backed MySQL emulation (embd_data.mysql_emulation; the
database file is ``$EMBD_SQLITE``), so the original SQL runs unchanged and writes stay local.
Values written into SQL text are formatted as Python 2 did (``py2_str``), so stored numbers are
identical to the legacy run. Functions on the per-embryo ``embryos`` table are out of scope.
"""
import json

import numpy as np

from embd_data import mysql_emulation as MySQLdb

from ._compat import not_ported, py2_str
from .varLookup import PARAM_NAMES

get_embryo_id = not_ported('db_utils_embryos.get_embryo_id')
add_float_column_2embryos = not_ported('db_utils_embryos.add_float_column_2embryos')
add_embryo = not_ported('db_utils_embryos.add_embryo')
get_tot_num_embs = not_ported('db_utils_embryos.get_tot_num_embs')
get_all_params = not_ported('db_utils_embryos.get_all_params')
get_embryo_number_by_rna = not_ported('db_utils_embryos.get_embryo_number_by_rna')
get_file_path_by_rna = not_ported('db_utils_embryos.get_file_path_by_rna')


def initialize_db():
    """
    initializes the database
    :return: connection and cursor to the database
    """
    myConn = MySQLdb.connect()
    return myConn, myConn.cursor()


def print_rows(cursor, limit=None):
    """
    prints all genes present in the database
    :return:
    """
    sql = "SELECT * FROM genes"
    if limit is not None:
        sql += " LIMIT {}".format(limit)
    cursor.execute(sql)
    rows = cursor.fetchall()
    for row in rows:
        print(row)


def get_columns(cursor, table):
    """
    print table column names
    :return:
    """
    cursor.execute("SHOW COLUMNS FROM {t}".format(t=table))
    result = [value[0] if value[0] != 'bY' else '`bY`' for value in cursor.fetchall()]
    return result


def get_all_ids(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT id FROM genes"
    cursor.execute(sql)
    return list(np.array(cursor.fetchall()).ravel())


def create_param_names(cursor):
    """
    adds parameters to the database
    :return:
    """

    columns = get_columns(cursor, 'genes')
    p_name_prev = ''
    for p_name in PARAM_NAMES:
        print(p_name)
        if p_name + '_GLS' not in columns:
            sql = "ALTER TABLE genes ADD COLUMN {} FLOAT NULL".format(p_name + '_GLS')
            if p_name_prev != '':
                sql += " AFTER {}".format(p_name_prev + '_GLS')
            cursor.execute(sql)
        if p_name + '_MS' not in columns:
            sql = "ALTER TABLE genes ADD COLUMN {} FLOAT NULL".format(p_name + '_MS')
            if p_name_prev != '':
                sql += " AFTER {}".format(p_name_prev + '_MS')
            cursor.execute(sql)


def insert_row(colms, cursor, table):
    """
    adds a new row to specified table in mySQL database
    :param colms: dictionary with column names and values
    :param cursor: cursor
    :param table: db table
    :return: None
    """

    names = ""
    values = ""
    my_colms = get_columns(cursor, table)
    for col in colms:  # cols is the keyname for key value pair in dict
        if col in my_colms and colms[col] is not None:  # colms[col] is the value for key value pair in dict
            if isinstance(colms[col], str):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=colms[col])
            elif isinstance(colms[col], np.ndarray):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=json.dumps(colms[col].tolist()))
            elif isinstance(colms[col], (list, dict)):
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=json.dumps(colms[col]))
            elif not np.isnan(colms[col]):  # if a single value (float or int)
                names += '{pn}, '.format(pn=col)
                values += '\'{pv}\', '.format(pv=py2_str(colms[col]))  # PORT: py2 str() rendering
            elif np.isnan(colms[col]):
                names += '{pn}, '.format(pn=col)
                values += 'NULL, '.format(pv=colms[col])
    names = names[:-2]  # removes the final comma and space from the string
    values = values[:-2]
    sql = "INSERT INTO {t} ({n}) VALUES ({v})".format(t=table, n=names, v=values)
    cursor.execute(sql)


def get_row(myid, cursor, table, colms=None, sonn=False, germ=False, ns=False):
    """
    reads columns in a row from mysql database and returns a dictionary
    :param myid: row id (int)
    :param cursor: mysql cursor
    :param table: table name (str)
    :param colms: list of names of columns
    :return: dictionary
    """

    sql = "SELECT "
    my_colms = get_columns(cursor, table)
    if colms is None:
        sql += '*  '
    else:
        for col in colms:
            if col in my_colms:
                sql += '{c}, '.format(c=col)
    if sonn:
        sql = sql[:-2] + ' FROM {t} WHERE id_sonn_key={id}'.format(t=table, id=myid)
    elif germ:
        sql = sql[:-2] + ' FROM {t} WHERE id_germ_key={id}'.format(t=table, id=myid)
    elif ns:
        sql = sql[:-2] + ' FROM {t} WHERE id_ns_key={id}'.format(t=table, id=myid)
    else:
        sql = sql[:-2] + ' FROM {t} WHERE id={id}'.format(t=table, id=myid)

    cursor.execute(sql)
    vals = list(cursor.fetchall()[0])
    for j in range(len(vals)):
        if isinstance(vals[j], str) and ('[' in vals[j] or '{' in vals[j]):
            vals[j] = json.loads(vals[j])
            if isinstance(vals[j], list):
                vals[j] = np.array(vals[j])
        elif vals[j] is None:
            vals[j] = np.nan
    field_names = [d[0] for d in cursor.description]
    return dict(zip(field_names, vals))


def update_row(myid, colms, cursor, table, sonn=False, germ=False):
    """
    update columns in a given row
    :param myid: row id
    :param colms: dictionary with column names and values
    :param cursor: cursor
    :param table: db table
    :return:
    """

    sql = "UPDATE {t} SET ".format(t=table)
    my_colms = get_columns(cursor, table)
    for col in colms:
        if col in my_colms and colms[col] is not None:
            if isinstance(colms[col], str):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=colms[col])
            elif isinstance(colms[col], np.ndarray):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=json.dumps(colms[col].tolist()))
            elif isinstance(colms[col], (list, dict)):
                sql += '{pn}=\'{pv}\', '.format(pn=col, pv=json.dumps(colms[col]))
            elif not np.isnan(colms[col]):
                sql += '{pn}={pv}, '.format(pn=col, pv=py2_str(colms[col]))  # PORT: py2 str() rendering
            elif np.isnan(colms[col]):
                sql += '{pn}=NULL, '.format(pn=col)
    if sonn:
        sql = sql[:-2] + ' WHERE id_sonn_key={id}'.format(id=myid)
    elif germ:
        sql = sql[:-2] + ' WHERE id_germ_key={id}'.format(id=myid)
    else:
        sql = sql[:-2] + ' WHERE id={id}'.format(id=myid)
    cursor.execute(sql)

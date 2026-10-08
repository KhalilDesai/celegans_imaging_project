"""Data layer for the EMBD downstream analysis.

Serves the published per-gene parameter table (``data/All_genes_mySQL_table.xlsx``, a dump of
the lab's MySQL ``genes`` table) in place of the MySQL database. It never writes to a database
server. Writes from the original code go to a local SQLite file.

* ``sheet``           - xlsx -> canonical CSV (exact float round-trip), gene_key reconstruction.
                         Python 3 only (needs openpyxl).
* ``mysql_emulation`` - a ``MySQLdb``-compatible ``connect()`` backed by SQLite. Runs unmodified
                         under Python 2.7 and Python 3.
"""

"""
Port of development/db_utils_genes.py: the functions on the downstream pipeline.

Changes (DECISIONS L4/L6):
* functions that used the ``__main__``-level global ``conn`` take ``conn`` as a parameter
  (default: the cursor's own connection);
* the hand-edited loop bounds of ``update_top_hits_PAD`` (``[1406, 674, ...]``) and
  ``update_column_reciprocal_matches`` (``range(310, 311)``) are an ``embd_ids`` parameter whose
  default is the original value;
* ``print('...').format(i)`` (a py2 print statement) is ``print('...'.format(i))``.
Utilities that populate MySQL from CSV files absent from the public data are not ported.
"""
import numpy as np
import pandas as pd

from . import db_utils_embryos


def _needs_lab_files(name):
    def _f(*args, **kwargs):
        raise NotImplementedError('%s: not ported (needs input files that are not in the public data)' % name)
    _f.__name__ = name
    return _f


for _n in ('print_rows_ce_id', 'check_gene_present', 'get_cancer_rank', 'update_cancer_rank',
           'update_human_ortho2_germ', 'update_embd_info', 'update_interacting_genes', 'read_manual_data',
           'update_manual_data', 'update_gene_wb_info', 'read_wb_info', 'read_embd_info', 'read_cancer_rank',
           'get_ensemble_ref', 'update_other_names', 'remove_gene', 'force_version', 'get_embryo_count_by_RNAi',
           'get_id_from_sonn', 'get_id_from_germ', 'get_id_from_ns'):
    globals()[_n] = _needs_lab_files('db_utils_genes.' + _n)


def _conn_of(cursor, conn):
    return conn if conn is not None else cursor._conn


def get_columns(cursor, table):
    """
    print table column names
    :return:
    """
    cursor.execute("SHOW COLUMNS FROM {t}".format(t=table))
    result = [value[0] for value in cursor.fetchall()]
    return result


def get_tot_num_genes(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT COUNT(*) FROM genes"
    cursor.execute(sql)
    return cursor.fetchall()[0][0]


def get_id_from_embd(embd_id, cursor):
    """
    returns id for first gene with embd number
    :param embd_id: embd number
    :param cursor: db cursor
    :return: int unique identifier for mySQL table
    """

    if isinstance(embd_id, str):  # PORT: py2 basestring
        embd_id = int(embd_id[-4:])
    sql = "SELECT id FROM genes WHERE embd_id={}".format(embd_id)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]


def get_id_from_genename(genename, cursor):
    """
    returns id for first gene with genename (i.e.'R05H5.2')
    :param embd_id: genename
    :param cursor: db cursor
    :return: int unique identifier for mySQL table
    """

    sql = "SELECT id FROM genes WHERE ce_id=\"{}\"".format(genename)
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    if len(sql_rows) == 0:
        return None
    else:
        return sql_rows[0][0]


def get_max_id(cursor):
    """
    calculates total number of gene entries in the database
    :return: int
    """

    sql = "SELECT id FROM genes ORDER BY id DESC LIMIT 0, 1"
    cursor.execute(sql)
    return cursor.fetchall()[0][0]


def create_columns_pad_rank(cursor):
    """
    adds columns to pad_rank database for top 50 hits (rank#) and pad values for the top 50 hits (pad#)
    - use this to build the pad_rank database
    :return:
    """

    columns = get_columns(cursor, 'pad_rank')  # retrieves columns from table as list
    for i in range(1, 51):  # all potential RNAi conditions in the EMBD screen
        print(i)
        if i not in columns:
            sql = "ALTER TABLE pad_rank ADD COLUMN `{}` INT(11) NULL".format('rank' + '{}'.format(i))  # adds columns for the top 50 RNAi
            # conditions (result formatted as INT)
            cursor.execute(sql)
            sql = "ALTER TABLE pad_rank ADD COLUMN `{}` FLOAT NULL".format('pad' + '{}'.format(i))# adds columns for PAD values for the
            # top 50 RNAi conditions (result formatted as float)
            cursor.execute(sql)


def get_reciprocal_matches(seed, top_hit_number, cursor, thr=0.55):
    """
    Gets top PAD hits from SQL and checks to see if the hit is reciprocal. Criteria is top x hits in first direction, top
    x hits in reverse direction... as long as PAD values are above 0.55 (or whatever thresh is established)
    :param seed: EMBD number
    :param top_hit_number: the number of top hits to consider in first direction (x)
    :param thresh- adjustable to determine what threshold of PAD to consider for best hits
    :return: best hit list (a list of RNAi conditions that show reciprocal top hits)
    """
    import numpy as np

    ranks_considered = []  # generates a list of 'rank1, rank2, rank3...' based on the specified top hit number
    pad_considered = []  # generates a list of 'pad1, pad2, pad3...' based on top hit number

    for i in range(1, top_hit_number+1):  # generates the two lists above
        rank = 'rank{}'.format(i)
        ranks_considered.append(rank)
        pad_val = 'pad{}'.format(i)
        pad_considered.append(pad_val)

    ranks_considered_string = ','.join(ranks_considered)  # converts ranks_considered to a string for SQL query
    pad_considered_string = ','.join(pad_considered)  # converts pad_considered to a string for SQL query


    sql = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(ranks_considered_string, seed)
    sql2 = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(pad_considered_string, seed)

    cursor.execute(sql)
    hits = cursor.fetchall()
    hits = np.array(hits[0])  # convert retrieved dataframe to np.array
    cursor.execute(sql2)
    pad = cursor.fetchall()
    pad = np.array(pad[0])

    high_pad = pad > thr  # generates a TRUE/FALSE list to ID usable indices
    inds = np.where(high_pad==True)  # retrieves inds where high_pad is true to apply to ranks considered
    ranks_use = hits[inds]  # a new list of ranks to use in the query that have PADs >0.5

    recip_hit_number = top_hit_number  # generates a reciprocal rank value (if top hit number is 10, this will be also be 10)
    # recip_hit_number = top_hit_number * 2  # generates a reciprocal rank value (if top hit number is 10, this will be 20)
    recip_ranks_considered = []  # generates a list of 'rank1, rank2, rank3...' based on the specified top hit number

    for j in range(1, recip_hit_number+1):  # generates the two lists above #This was missing the +1 in the first iteration- is fixed
        recip_rank = 'rank{}'.format(j)
        recip_ranks_considered.append(recip_rank)
    recip_ranks_considered_string = ','.join(recip_ranks_considered)

    best_hits = []  # list of hits that are reciprocal hits. Populated by loop below
    for hit in ranks_use:
        sql3 = "SELECT {0} FROM PAD_RANK WHERE rnai={1};".format(recip_ranks_considered_string, hit)
        cursor.execute(sql3)
        recip_hits = cursor.fetchall()
        recip_hits = recip_hits[0]
        if seed in recip_hits:
            best_hits.append(hit)
    print((seed, best_hits))
    return(best_hits)


def update_automated_description():
    '''
    Populates automated description column based on proximity of gene of interest to centroid of manual groups
    :param cursor:
    :return:
    '''
    from .groupAssesment import get_all_closest_groups

    get_all_closest_groups(to_csv=False, to_SQL=True, to_plot=False)


ORIGINAL_TOP_HITS_LIST = [1406, 674, 952, 866, 885, 703, 1360, 1904, 1903, 1905]


def update_top_hits_PAD(cursor, conn=None, embd_ids=ORIGINAL_TOP_HITS_LIST):
    '''
    Populates top_PAD_hits column in mySQL with top 10 PAD hits. All hits must have a PAD above 0.5.
    :param cursor:
    :return:
    '''
    import numpy as np
    from .embdFunc import getGeneCommonName
    conn = _conn_of(cursor, conn)  # PORT: the original used the __main__-level global `conn`
    thr = 0.5

    # for i in range(1,504):
    for i in embd_ids:  # PORT: the hand-edited list [1406,674,952,866,885,703,1360,1904,1903,1905] is the default

        print('starting {0}'.format(i))  # PORT: py2 print statement
        seed = i
        ranks_considered = 'rank1, rank2, rank3, rank4, rank5, rank6, rank7, rank8, rank9, rank10'
        pad_considered = 'pad1, pad2, pad3, pad4, pad5, pad6, pad7, pad8, pad9, pad10'
        sql = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(ranks_considered, seed)
        sql2 = 'SELECT {0} FROM PAD_RANK WHERE rnai={1};'.format(pad_considered, seed)

        cursor.execute(sql)  # gets EMBD top 10
        hits = cursor.fetchall()
        hits = np.array(hits[0])  # convert retrieved dataframe to np.array

        cursor.execute(sql2)  # gets PADs, top 10
        pad = cursor.fetchall()
        pad = np.array(pad[0])

        high_pad = pad > thr  # generates a TRUE/FALSE list to ID usable indices (above 0.5 PAD)
        inds = np.where(high_pad==True)  # retrieves inds where high_pad is true to apply to ranks considered
        ranks_use = hits[inds]  # a new list of ranks to use in the query that have PADs >thresh

        genes = []  # populated by common gene name for each EMBD_id
        for item in ranks_use:
            genes.append(getGeneCommonName(item))
        genes_use = ', '.join(genes)  # makes a string to insert into mySQL query

        sql = "SELECT id, embd_id, top_hits_PAD FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, top_hits_PAD = cursor.fetchall()[0]
        sql2 = "UPDATE genes SET top_hits_PAD=\"{h}\" WHERE id={id};".format(h=genes_use, id=id)
        cursor.execute(sql2)
        conn.commit()


def update_WT_like_or_phenotype(cursor, conn=None):
    '''
    Populates column in SQL to define if RNAi condition is called as 'WT-like' (moving and sustaining movement), or
    'phenotype' (movement is 0 for either GLS, MS or both). Reads movementGLS and movementMS to determine.
    :param cursor:
    :return: populates mySQL
    '''
    conn = _conn_of(cursor, conn)  # PORT: the original used the __main__-level global `conn`
    for i in range(1,504):
        print('starting {0}'.format(i))  # PORT: py2 print statement
        sql = "SELECT id, embd_id, movement_GLS, movement_MS, auto_description, WT_like_or_phenotype FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, movementG, movementM, auto_description, WT_like_or_phenotype = cursor.fetchall()[0]
        if movementG == 0 or movementM == 0:
            sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='Phenotype', id=id)
            cursor.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1:
            if "Embryos develop at a slow pace." in auto_description:
                sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='Dev. Delay', id=id)
                cursor.execute(sql2)
                conn.commit()
            else:
                sql2 = "UPDATE genes SET WT_like_or_phenotype=\"{p}\" WHERE id={id};".format(p='WT-like', id=id)
                cursor.execute(sql2)
                conn.commit()


def update_column_reciprocal_matches(cursor, conn=None, embd_ids=range(310, 311)):
    """
    populates reciprocal matches column in genes database
    :return:
    """
    from .embdFunc import getGeneCommonName
    conn = _conn_of(cursor, conn)  # PORT: the original used the __main__-level global `conn`

    for i in embd_ids:  # PORT: the hand-edited range(310,311) is the default
        best_hits = get_reciprocal_matches(i, 10, cursor, thr=0.55)
        top_five = best_hits[0:10]  # changed to top 10

        genes = []
        for item in best_hits:
            genes.append(getGeneCommonName(item))
        top_five_genes = ', '.join(genes[0:10])
        print((i,top_five, top_five_genes))

        sql = "SELECT id, embd_id, movement_GLS, movement_MS, WT_like_or_phenotype, reciprocal_matches FROM genes WHERE embd_id={}".format(i)
        cursor.execute(sql)
        id, embd_id, movementG, movementM, WT_like_or_phenotype, reciprocal_matches = cursor.fetchall()[0]

        if movementG == 0 or movementM == 0:  # checks to see if there is a scored phenotype by checking movement
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        elif movementG == 1 and movementM == 1 and "Dev. Delay" in WT_like_or_phenotype:  # looks for dev. delay
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=top_five_genes, id=id)
            cursor.execute(sql2)
            conn.commit()
        else:  # if no scored phenotype (movement) or dev delay, doesnt provide 'hits' as it is expected to be wild-type.
            wt_note = 'Embryos do not exhibit an arrest phenotype'
            sql2 = "UPDATE genes SET reciprocal_matches=\"{bh}\" WHERE id={id};".format(bh=wt_note, id=id)
            cursor.execute(sql2)
            conn.commit()


def get_sql_data_for_man_group_conditions():
    '''
    retrieves the automated parameters and manual parameters for the conditions selected in manual groupings from the SQL. To be used for classifier algorithms.
    :param cursor:
    :return: df of automated and manual params for 94 manually grouped conditions
    '''
    import numpy as np
    conn, cursor = db_utils_embryos.initialize_db()
    sql = "SELECT id, embd_id, man_group, movement_GLS, movement_MS, WT_like_or_phenotype, maxG_GLS, aG_GLS, sY_GLS, fracR_GLS, fracG_GLS, CoM0Rtd_GLS_GLS, CoM1Rscale_GLS_GLS, MoI0GavgRes_GLS_GLS, MoI1Gscale_GLS_GLS, " \
          "MoI1GstdTail_GLS_GLS, MoI0Rtd_GLS_GLS, MoI0RavgRes_GLS_GLS, MoI1Rtd_GLS_GLS, MoI0Ytd_GLS_GLS, MoI0YavgRes_GLS_GLS, MoI1Yscale_GLS_GLS, MoI1YavgRes_GLS_GLS, sG_MS, maxR_MS, aR_MS, sR_MS, maxHead_MS, tailHead_MS, " \
          "devHead_MS, devLength_MS, tailLength_MS, CoM0Gtd_MS_MS, CoM0Gscale_MS_MS, CoM1Gscale_MS_MS, CoM1GavgRes_MS_MS, CoM0Rscale_MS_MS, CoM0RavgRes_MS_MS, CoM1Rtd_MS_MS, CoM1RavgRes_MS_MS, CoM1RstdTail_MS_MS, " \
          "MoI0Gscale_MS_MS, MoI1Gscale_MS_MS, MoI1GstdTail_MS_MS, MoI1Rscale_MS_MS, MoI1RavgRes_MS_MS, MoI1RstdTail_MS_MS,man_GS1, man_GS2, man_GS3, man_GS4, man_GS5, man_GS6,man_GS7, man_GS8, man_GS9, man_GS12, man_GS14," \
          "man_GS15, man_GS19, man_GS20, man_GS21, man_GS22, man_GSemb_comments, man_MS1, man_MS2, man_MS4, man_MS5, man_MS6, man_MS7, man_MS8, man_MS9, man_MS10, man_MS11, man_MS13, man_MS14, man_MS15, man_MS17, man_MS19," \
          "man_MS20, man_MS21, man_MS22, man_MSemb_comments FROM genes WHERE man_group IS NOT NULL"
    cursor.execute(sql)
    df = pd.DataFrame(cursor.fetchall(), columns=[d[0] for d in cursor.description])  # PORT: was pd.read_sql_query(sql, conn), which needs a real DB-API/SQLAlchemy connection
     #note that fetchall doesnt retrieve a df and this command works better
    return df

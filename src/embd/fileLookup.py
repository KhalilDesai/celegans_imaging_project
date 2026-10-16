"""Port of development/fileLookup.py: same names, rooted at config.folder_in() instead of 'Z:/'."""
from .varLookup import FOLDER_IN

FILENAME_TRACKING_SHEET = FOLDER_IN + 'EMBD_fileNames_Tracking_Sheet.csv'
FILENAME_PAD = FOLDER_IN + 'phenoDist_rank_04142022.csv'  # Includes fixed arx-1 and arx-3 and fixed headInt params as of 4/14/22

FILENAME_NETWORK_DATA = FOLDER_IN + 'phenoDistCalc/phenoDist.csv'
FOLDER_RNAI_DATA = FOLDER_IN + 'RNAi_Data/'
FOLDER_TIME_ALIGN_MOVIES = FOLDER_IN + 'Time_align_movies/'
FOLDER_TOP_HITS_MOVIES = FOLDER_IN + 'Top_hits/'

FILENAME_GENE_NAMES = FOLDER_IN + 'gene_key.csv'
FILENAME_COMPLETE_GENE_LIST = FOLDER_IN + 'EMBD_fullGenelist.csv'
FILENAME_GENE_WB_INFO = FOLDER_IN + 'genes_wb_info3b.csv'
FILENAME_EMBD_INFO = FOLDER_IN + 'EMBD_gene_info.csv'
FILENAME_GENE_CANCER_RANK = FOLDER_IN + 'Divoli_TS_OG_ranks.csv'
FILENAME_MANUAL_SCORING = FOLDER_IN + 'manualScoringForMySQL_2020update.csv'
FILENAME_EMBD_ORTHOLIST = FOLDER_IN + 'ortholist_results_embd.csv'
FILENAME_OTHER_NAMES = FOLDER_IN + 'other_names_list_080422.csv'

FILENAME_WEB_CONTENT = FOLDER_IN + 'WebContent/'

FILENAME_SONNICHSEN = FOLDER_IN + 'prior_screen_data/sonnichsen_gene_info_table.csv'
FILENAME_SONNICHSEN_CSI = FOLDER_IN + 'prior_screen_data/Sonnichsen_data_CSI.csv'
FILENAME_SONNICHSEN_ORTHOLIST = FOLDER_IN + 'Gene_lists/ortholist_results_sonnichsen.csv'

FILENAME_GERMLINE = FOLDER_IN + 'prior_screen_data/green2011_gene_info_table3.csv'
FILENAME_GERMLINE_CSI = FOLDER_IN + 'prior_screen_data/Green2011_data_CSI.csv'
FILENAME_GERMLINE_ORTHOLIST = FOLDER_IN + 'Gene_lists/ortholist_results_germline.csv'

# Not defined in the public repository but imported by four functions (DECISIONS C5, G6, Q4).
# Empty until the lab provides the original list.
SKIPPED = []

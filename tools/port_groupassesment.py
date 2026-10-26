"""Generate src/embd/groupAssesment.py from development/groupAssesment.py (published-path functions).

Copies the original function bodies and applies only the listed edits, so the port can be
regenerated and audited. Run: python tools/port_groupassesment.py
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KEEP = ['getDist2WT', 'getDist2EO', 'getDist2AVG', 'getAvgVector', 'getPAD2EO', 'get_all_PAD2EO', 'get_all_dist2EO',
        'calculate_avg_vectors_man_groups', 'get_manual_groups', 'get_rnai_groups_pheno_distinct',
        'find_closest_groups', 'get_closest_group_description', 'get_group_description',
        'get_all_closest_groups', 'make_UMAP_content_files']
EDITS = [
    # NOTE (comment only, no behaviour change; DECISIONS A17): duplicated genes in the manual groups.
    ("    stable list to access the 17 groups genes\n    :return: rnaiList, labels, defect_label\n    \'\'\'\n",
     "    stable list to access the 17 groups genes\n    :return: rnaiList, labels, defect_label\n    \'\'\'\n"
     "    # NOTE (port, DECISIONS A17): three genes were knocked down twice in EMBD 1-503, most likely by accident\n"
     "    # (it was not clear at the time that they were the same gene): elt-1 = EMBD 16 & 34, vab-1 = 27 & 73,\n"
     "    # apx-1 = 31 & 315. The original code never picks one: both are separate conditions everywhere. Here the\n"
     "    # 'dorsRupt' group [34, 16, 31, 264] contains BOTH elt-1 conditions (two replicates of one gene), which\n"
     "    # makes that group look tighter and slightly biases the group-based parameter optimisation and the\n"
     "    # centroid descriptions. 315, 27 and 73 are in no group. Kept as published; not changed pending a decision.\n"),
    # PORT (perf, PERFORMANCE #4): all-pairs PAD / distance rows computed by _fast.PairwiseRows, which is
    # bit-identical to the per-pair RNAiClass.getPAD / getDistance calls it replaces.
    ("    pad = []\n    for i in range(len(rna_list)):\n        for j in range(len(rna_list)):\n            p = rna_list[i].getPAD(rna_list[j])",
     "    pad = []\n    from ._fast import PairwiseRows  # PORT (perf #4)\n    _rows = PairwiseRows(rna_list)  # PORT (perf #4)\n"
     "    for i in range(len(rna_list)):\n        _row_i = _rows.row(i)[1]  # PORT (perf #4): == [rna_list[i].getPAD(r) for r in rna_list], bit-identical\n"
     "        for j in range(len(rna_list)):\n            p = _row_i[j]  # original: rna_list[i].getPAD(rna_list[j])"),
    ("    dist = []\n    for i in range(len(rna_list)):\n        for j in range(len(rna_list)):\n            p = rna_list[i].getDistance(\n                rna_list[j])",
     "    dist = []\n    from ._fast import PairwiseRows  # PORT (perf #4)\n    _rows = PairwiseRows(rna_list)  # PORT (perf #4)\n"
     "    for i in range(len(rna_list)):\n        _row_i = _rows.row(i)[0]  # PORT (perf #4): == [rna_list[i].getDistance(r) for r in rna_list], bit-identical\n"
     "        for j in range(len(rna_list)):\n            p = _row_i[j]  # original: rna_list[i].getDistance(rna_list[j])"),
    ("    print np.array(pad)\n", "    print(np.array(pad))\n"),
    ("        if pad > 0:\n            print pad\n", "        if pad > 0:\n            print(pad)\n"),
    ("group [0] is not 'dev delay'", "group [0] != 'dev delay'"),
    ("    from embdFunc import getGeneCommonName, getNDDistance\n", "    from .embdFunc import getGeneCommonName, getNDDistance\n"),
    ("        from db_utils_embryos import initialize_db\n", "        from .db_utils_embryos import initialize_db\n"),
    ("    from fileLookup import SKIPPED\n", "    from .fileLookup import SKIPPED\n"),
    ("    path = 'Z:/closest_groups{}.csv'.format(time.strftime('%Y%m%d'))",
     "    path = FOLDER_IN + 'closest_groups{}.csv'.format(time.strftime('%Y%m%d'))"),
    ("    np.save('Z:/PAD_Matrix_msParams.npy', pad_array_reshaped)\n    np.save('Z:/Dist_Matrix_msParams.npy', dist_array_reshaped)\n    np.save('Z:/Matrix_labels_msParams.npy', labels)",
     "    np.save(FOLDER_IN + 'PAD_Matrix_msParams.npy', pad_array_reshaped)\n    np.save(FOLDER_IN + 'Dist_Matrix_msParams.npy', dist_array_reshaped)\n    np.save(FOLDER_IN + 'Matrix_labels_msParams.npy', labels)"),
    ("    np.save('Z:/group_Lists_msParams.npy', groupLists)\n",
     "    np.save(FOLDER_IN + 'group_Lists_msParams.npy', np.array(groupLists, dtype=object), allow_pickle=True)  # ragged (B12)\n"),
    ("    np.save('Z:/group_Labels_msParams.npy', groupLabels)\n    np.save('Z:/group_colors_msParams.npy', groupColors)\n    np.save('Z:/group_Labels_toPlot_msParams.npy', plot_labels)",
     "    np.save(FOLDER_IN + 'group_Labels_msParams.npy', groupLabels)\n    np.save(FOLDER_IN + 'group_colors_msParams.npy', groupColors)\n    np.save(FOLDER_IN + 'group_Labels_toPlot_msParams.npy', plot_labels)"),
    ("    from sklearn.manifold import TSNE\n    import pandas as pd\n    import time\n    from embdFunc import getGeneCommonName\n\n    rnaiList",
     "    rnaiList"),
]
HEADER = '"""\nPort of development/groupAssesment.py (Created on Mar 3, 2017): the functions on the\npublished-output path: manual groups and centroids, closest group / auto_description,\nall-pairs PAD and distance matrices, UMAP input files.\n\nGENERATED by tools/port_groupassesment.py, which copies the original function bodies and applies\nonly these edits: py2 print statements -> print(); ``is not \'dev delay\'`` -> ``!=`` (same result\nhere); package-relative imports; \'Z:/\' -> FOLDER_IN; the ragged group list is saved as an\nexplicit object array (numpy >= 1.24, B12); unused function-local imports in\nmake_UMAP_content_files dropped.\n\nNOT YET PORTED (not on the path to any published output; analysis/optimisation experiments,\nseveral already broken in the original, DECISIONS G7): getCSI2EO, plot*, getRNAiGroupsArrest,\ngetRNAiGroupsPhenoOld, getRNAiGroupsPheno, retrieve_and_map_man_group_autodata,\ntest_two_sample_ks_test, log_regression, remove_param_from_use, optimizeByGain,\nget_rnai_groups_obj, optimize_groups_by_min_max, optimizeByInteractome, flatten_groups,\noptimize_groups_by_rank, calc_pad_normdist, get_top_hits, get_in_group_rank, evaluate_ranking,\ncross_corr_params, get_PAD_matrix, get_dist_matrix, get_tSNE, show_tsne, verify_tsne,\nshow_clustermap, show_heatmap, check_dist2zero, find_closest_group_byEMBD,\nplot_dist2WT_dist2avg, make_tSNE_all_data.\n"""\nimport numpy as np\n\nfrom . import embdFunc\nfrom .params_use_assign import get_params_use\nfrom .RNAiClass import RNAiClass\nfrom .varLookup import FOLDER_IN\n\nPARAM_NAMES_USE_GLS, PARAM_NAMES_USE_MS = get_params_use()\n\n\n'


def body(src, name):
    idx = src.index('\ndef %s(' % name)
    end = src.find('\ndef ', idx + 1)
    return src[idx + 1: end if end != -1 else len(src)].rstrip() + '\n'


def main():
    src = open(os.path.join(REPO, 'development', 'groupAssesment.py'), encoding='latin-1').read()
    code = '\n\n'.join(body(src, n) for n in KEEP)
    for a, b in EDITS:
        assert code.count(a) == 1, 'edit anchor not found exactly once: %r' % a[:60]
        code = code.replace(a, b)
    for line in code.splitlines():
        s = line.split('#', 1)[0]
        assert not s.lstrip().startswith('print ') and ' print ' not in s, 'py2 print left: %r' % line
    with open(os.path.join(REPO, 'src', 'embd', 'groupAssesment.py'), 'w') as f:
        f.write(HEADER + code)


if __name__ == '__main__':
    main()

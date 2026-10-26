"""
Port of development/makeUMAPs_V2.py (originally a module-level script run under Python 3.6).

Wrapped in ``run()``: ``np.random.seed(71076)`` then ``umap.UMAP()`` with the original defaults, the
503x503 PAD matrix rows as features, the
11 manual groups + query overlay, and the same labels/colours/offsets. The input matrix and output
paths are parameters (the original hardcoded 'Z://PAD_Matrix_042122_optParams.npy' and an SVG
name). Returns the embedding.

Reproducibility (maintainer decision 2026-10-26, NUMERICAL_DIFFERENCES N3): by default UMAP is seeded
with ``random_state=71076`` (umap-learn then runs single-threaded), so repeated runs give identical
embeddings. ``random_state=None`` reproduces the original call (unseeded, parallel, not reproducible).
"""
import numpy as np


def run(data, query=1, common_name=None, svg_path=None, title='UMAP of PAD with optimal parameters',
        seed=71076, random_state=71076):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import umap

    sns.set(style='white', context='poster', rc={'figure.figsize': (14, 10)})
    np.random.seed(seed)

    rnaiList = [i for i in range(1, 504)]
    labels = [str(x) for x in rnaiList]

    groupLists = [[21, 364, 408, 130], [417, 64, 115, 281, 77], [77, 63, 19, 501],
                  [264, 386, 31, 32, 255, 388, 422, 118, 359],
                  [181, 182, 357, 184, 185, 154, 363, 117, 45, 447, 108], [9, 398, 435, 45, 52, 38],
                  [95, 4, 57, 5, 277],
                  [26, 453, 235, 327, 489, 379, 420, 225, 289, 261, 350],
                  [10, 15, 217, 18, 28, 439, 177, 291, 209], [110, 142, 98, 101, 186],
                  [383, 495, 498, 414, 375, 396, 321]]
    groupLists += [[query]]
    groupColors = ['sienna', 'orangered', 'firebrick', 'darkorange', 'gold', 'y', 'green', 'turquoise', 'teal',
                   'dodgerblue', 'b', 'deeppink']
    plot_labels = ['Early Arrest, \n High Ectoderm', 'Early Arrest, \n High Mesoderm', 'Mid Arrest, \n High Mesoderm',
                   'Early Arrest Sectored  \n Germ-Layers', 'Elongation Defect,  \n High Neurons',
                   'Enclosure/Early Elong.  \n Dorsal Defect', 'Enclosure/Early Elong. \n Ventral Rupture',
                   'Variable \n Defects', 'Late Elongation \n Arrest', 'Dev. Delay', 'Wild-type-like', 'Query']

    fit = umap.UMAP() if random_state is None else umap.UMAP(random_state=random_state)  # original: umap.UMAP()
    points = fit.fit_transform(data)

    if svg_path is not None:
        plt.figure(figsize=(12, 12))
        for point in points:
            plt.scatter(point[0], point[1], color="grey", marker=".", alpha=0.5)
        for g in range(len(groupLists)):
            man_group_list = groupLists[g]
            lab = plot_labels[g]
            inds1 = [labels.index(str(p)) for p in man_group_list]
            for ind in inds1:
                plt.scatter(points[ind][0], points[ind][1], color=groupColors[g], marker=".", alpha=0.5)
                plt.scatter(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1], axis=0)[1], label=lab,
                            marker="*", alpha=0.5, color=groupColors[g])
                if lab == 'Query':
                    query_label = '{n}(RNAi)'.format(n=common_name)
                    plt.text(np.mean(points[inds1] + 0.5, axis=0)[0], np.mean(points[inds1], axis=0)[1], query_label,
                             color=groupColors[g], size=14, style='italic', weight='bold')
                elif lab == 'Late Elongation \n Arrest':
                    plt.text(np.mean(points[inds1], axis=0)[0], np.mean(points[inds1] - 0.6, axis=0)[1], lab,
                             color=groupColors[g], size=11)
                elif lab == 'Mid Arrest, \n High Mesoderm':
                    plt.text(np.mean(points[inds1], axis=0)[0] + 0.1, np.mean(points[inds1] + 0.3, axis=0)[1], lab,
                             color=groupColors[g], size=11)
                elif (g % 2) == 0:
                    plt.text(np.mean(points[inds1] - 0.2, axis=0)[0], np.mean(points[inds1] - 0.5, axis=0)[1], lab,
                             color=groupColors[g], size=11)
                else:
                    plt.text(np.mean(points[inds1] + 0.2, axis=0)[0], np.mean(points[inds1] + 0.1, axis=0)[1], lab,
                             color=groupColors[g], size=11)
        plt.title(title)
        plt.savefig(svg_path)
        plt.close()
    return points

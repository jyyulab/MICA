#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from functools import partial
from sklearn.manifold import TSNE
import umap
from .distance import numba_calc_mi_dis


def visual_embed(clustering_res, frame_dr, out_dir, dis_metric='euclidean',
                 visual_method='UMAP', suffix=None, num_works=10, perplexity=30, min_dist=0.6):
    """ Visualize clustering results using tSNE for clustering results.
    Args:
        partition (dict): clustering results, {index: cluster label}
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
        dis_metric (str): distance metric (default: mi)
        visual_method (str): t-SNE or UMAP (default: UMAP)
        suffix (type): suffix to be added to output file name (default: None)
        perplexity (int): perplexity parameter in tSNE, Larger datasets usually require a larger
                          perplexity [default: 30]
        min_dist (float): min_dist parameter in UMAP, minimum distance of points in the embedded space
        num_works (int): number of threads for using numexpr package
    Outputs:
        TXT file with 2D embedded coordinates
        PNG image of 2D embedding
    """
    os.environ["NUMEXPR_MAX_THREADS"] = str(num_works)
    perplexity = np.min([perplexity, np.max(clustering_res.groupby(["label"]).size())])
    if dis_metric == 'mi':
        num_bins = int((frame_dr.shape[0]) ** (1 / 3.0))
        num_genes = frame_dr.shape[1]
        partial_calc_norm_mi = partial(numba_calc_mi_dis, bins=num_bins, m=num_genes)

        if visual_method == 't-SNE':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         metric=partial_calc_norm_mi, early_exaggeration=12.0, n_jobs=num_works).fit_transform(frame_dr)
        elif visual_method == 'UMAP':
            res = umap.UMAP(random_state=42, metric=partial_calc_norm_mi, n_neighbors=20,
                            min_dist=min_dist).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=clustering_res.index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid visualization method: {}'.format(visual_method))
    else:
        if visual_method == 't-SNE':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         early_exaggeration=12.0, n_jobs=num_works).fit_transform(frame_dr)
        elif visual_method == 'UMAP':
            res = umap.UMAP(random_state=42, metric='euclidean', n_neighbors=30,
                            min_dist=min_dist).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=clustering_res.index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid visualization method: {}'.format(visual_method))

    embed_2d = pd.DataFrame(data=embed, index=clustering_res.index, columns=["X", "Y"])
    clustering_res = pd.concat([clustering_res, embed_2d.loc[clustering_res.index, :]], axis=1)
    final_res = clustering_res.loc[:, ["X", "Y", "label"]]
    # save 2D embedding to txt file
    if suffix:
        out_txt_file = '{}/clustering_{}_{}_{}_{}.txt'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1],
                                                              suffix)
        out_png_file = '{}/clustering_{}_{}_{}_{}.png'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1],
                                                              suffix)
    else:
        out_txt_file = '{}/clustering_{}_{}_{}.txt'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1])
        out_png_file = '{}/clustering_{}_{}_{}.png'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1])
    final_res.to_csv(out_txt_file, sep="\t")
    scatter_plot(final_res, out_png_file, method=visual_method)


def scatter_plot(data, out_file, marker_size=1.0, marker="o", method='UMAP', marker_scale=10.0):
    if data is None:
        return
    plt.figure(figsize=(10, 10), dpi=400)
    lab = np.unique(data.loc[:, "label"])
    colors = plt.cm.jet(np.linspace(0, 1, len(lab)))
    for i, z in enumerate(lab):
        df = data.loc[data.loc[:, "label"] == z, ["X", "Y"]]
        plt.scatter(df.loc[:, "X"], df.loc[:, "Y"], facecolor=colors[i], s=marker_size,
                    marker=marker, vmin=0, vmax=len(lab), label=str(z) + "(" + str(df.shape[0]) + ")", alpha=0.7)
        center = np.mean(df, axis=0)
        plt.scatter(center.loc["X"], center.loc["Y"], marker="o", c="white", alpha=0.7, s=100, edgecolor="k")
        plt.scatter(center.loc["X"], center.loc["Y"], marker="$%d$" % (i+1), c="black", alpha=0.7, s=80, edgecolor="k")
    plt.ylabel("{}-2".format(method.upper()))
    plt.xlabel("{}-1".format(method.upper()))
    plt.xticks([])
    plt.yticks([])
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Clusters(" + str(data.shape[0]) + ")",
               markerscale=marker_scale)
    plt.savefig(out_file, bbox_inches="tight")
    plt.cla()

#!/usr/bin/env python3

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


def visual_embed_louvain(partition, resolution, index, frame_dr, out_dir, embed_method='umap', dis_metric='euclidean',
                         perplexity=30):
    """ Visualize clustering results using tSNE.
    Args:
        partition (dict): clustering results, {index: cluster label}
        resolution (float): Determines size of the communities
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
        embed_method (str): tsne or umap [default: umap]
        dis_metric (str): distance metric [default: mi]
        perplexity (int): perplexity parameter in tSNE, Larger datasets usually require a larger perplexity
                          [default: 30]
    Outputs:
        TXT file with 2D embedded coordinates
        PNG image of 2D embedding
    """
    labels = [x + 1 for x in partition.values()]
    clustering_res = pd.DataFrame(data=labels, index=index, columns=["label"])
    perplexity = np.min([perplexity, np.max(clustering_res.groupby(["label"]).size())])
    if dis_metric == 'mi':
        num_bins = int((frame_dr.shape[0]) ** (1 / 3.0))
        num_genes = frame_dr.shape[1]
        partial_calc_norm_mi = partial(numba_calc_mi_dis, bins=num_bins, m=num_genes)
        if embed_method == 'tsne':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         metric=partial_calc_norm_mi, early_exaggeration=12.0).fit_transform(frame_dr)
        elif embed_method == 'umap':
            res = umap.UMAP(random_state=42, metric=partial_calc_norm_mi, n_neighbors=20,
                            min_dist=0.2).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid embed method: {}'.format(embed_method))
    else:
        if embed_method == 'tsne':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         early_exaggeration=12.0).fit_transform(frame_dr)
        elif embed_method == 'umap':
            res = umap.UMAP(random_state=42, metric='euclidean', n_neighbors=20, min_dist=0.2).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid embed method: {}'.format(embed_method))
    embed_2d = pd.DataFrame(data=embed, index=index, columns=["X", "Y"])

    clustering_res = pd.concat([clustering_res, embed_2d.loc[clustering_res.index, :]], axis=1)
    final_res = clustering_res.loc[:, ["X", "Y", "label"]]
    # save 2D embedding to txt file
    out_txt_file = '{}/clustering_{}_{}_{}_{}.txt'.format(out_dir, embed_method, dis_metric, frame_dr.shape[1],
                                                          resolution)
    final_res.to_csv(out_txt_file, sep="\t")
    out_png_file = '{}/clustering_{}_{}_{}_{}.png'.format(out_dir, embed_method, dis_metric, frame_dr.shape[1],
                                                          resolution)
    scatter_plot(final_res, out_png_file)


def scatter_plot(data, out_file, marker_size=20, marker="o"):
    if data is None:
        return
    fig = plt.figure(figsize=(10, 10), dpi=300)
    lab = np.unique(data.loc[:, "label"])
    colors = plt.cm.jet(np.linspace(0, 1, len(lab)))
    for z in lab:
        df = data.loc[data.loc[:, "label"] == z, ["X", "Y"]]
        plt.scatter(df.loc[:, "X"], df.loc[:, "Y"], facecolor=colors[z-1], s=marker_size,
                    marker=marker, vmin=0, vmax=len(lab), label=str(z) + "(" + str(df.shape[0]) + ")", alpha=0.7)
        center = np.mean(df, axis=0)
        plt.scatter(center.loc["X"], center.loc["Y"], marker="o", c="white", alpha=0.7, s=100, edgecolor="k")
        plt.scatter(center.loc["X"], center.loc["Y"], marker="$%d$" % z, c="black", alpha=0.7, s=80, edgecolor="k")
    plt.ylabel("MICA-2")
    plt.xlabel("MICA-1")
    plt.xticks([])
    plt.yticks([])
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title="Clusters(" + str(data.shape[0]) + ")")
    plt.savefig(out_file, bbox_inches="tight")

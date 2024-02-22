#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from functools import partial
from sklearn.manifold import TSNE
import umap
from .distance import numba_calc_mi_dis
from sklearn.metrics import silhouette_samples


def visual_embed(clustering_res, frame_dr, out_dir, dr_dim=None, dis_metric='euclidean',
                 visual_method='UMAP', suffix=None, num_works=10, perplexity=30, min_dist=0.6, marker_size=1.0,
                 marker_scale=10.0):
    """ Visualize clustering results using tSNE for clustering results.
    Args:
        clustering_res (dict): clustering results, {sample index: cluster label}
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
        dr_dim (int): number of dimensions used for UMAP or tSNE embedding
        dis_metric (str): distance metric (default: mi)
        visual_method (str): t-SNE or UMAP (default: UMAP)
        suffix (type): suffix to be added to output file name (default: None)
        perplexity (int): perplexity parameter in tSNE, Larger datasets usually require a larger
                          perplexity [default: 30]
        min_dist (float): min_dist parameter in UMAP, minimum distance of points in the embedded space
        num_works (int): number of threads for using numexpr package
        marker_size (float): dot size in scatter plot
        marker_scale (float): dot size scale
    Outputs:
        TXT file with 2D embedded coordinates
        PNG image of 2D embedding
    """
    os.environ["NUMEXPR_MAX_THREADS"] = str(num_works)
    perplexity = np.min([perplexity, np.max(clustering_res.groupby(["label"]).size())])
    if dis_metric == 'mi':
        num_bins = int((frame_dr.shape[1]) ** (1 / 3.0))
        num_genes = frame_dr.shape[1]
        partial_calc_norm_mi = partial(numba_calc_mi_dis, bins=num_bins, m=num_genes)

        if visual_method.lower() == 't-sne':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         metric=partial_calc_norm_mi, early_exaggeration=12.0, n_jobs=num_works).fit_transform(frame_dr)
        elif visual_method.lower() == 'umap':
            res = umap.UMAP(random_state=42, metric=partial_calc_norm_mi, n_neighbors=20,
                            min_dist=min_dist).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=clustering_res.index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid visualization method: {}'.format(visual_method))
    else:
        if dr_dim:
            frame_dr = frame_dr.iloc[:, 0:dr_dim]
        if visual_method.lower() == 't-sne':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         early_exaggeration=12.0, n_jobs=num_works).fit_transform(frame_dr)
        elif visual_method.lower() == 'umap':
            res = umap.UMAP(random_state=42, metric='euclidean', n_neighbors=30,
                            min_dist=min_dist).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=clustering_res.index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid visualization method: {}'.format(visual_method))

    embed_2d = pd.DataFrame(data=embed, index=clustering_res.index, columns=["X", "Y"])
    clustering_res = pd.concat([clustering_res, embed_2d.loc[clustering_res.index, :]], axis=1)
    final_res = clustering_res.loc[:, ["X", "Y", "label"]]
    final_res.index.name = 'ID'
    # save 2D embedding to txt file
    if suffix:
        out_txt_file = '{}/clustering_{}_{}_{}_{}.txt'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1],
                                                              suffix)
        out_img_file = '{}/clustering_{}_{}_{}_{}.pdf'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1],
                                                              suffix)
    else:
        out_txt_file = '{}/clustering_{}_{}_{}.txt'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1])
        out_img_file = '{}/clustering_{}_{}_{}.pdf'.format(out_dir, visual_method, dis_metric, frame_dr.shape[1])
    final_res.to_csv(out_txt_file, sep="\t")
    scatter_plot(final_res, out_img_file, method=visual_method, marker_size=marker_size, marker_scale=marker_scale)
    return final_res


def scatter_plot(data, out_file, marker_size=1.0, marker="o", method='UMAP', marker_scale=10.0):
    """ Create a scatter plot for cluster results embedded in 2D space.
    """
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


def silhouette_plot(labels, frame_dr, num_clusters, silhouette_avg, out_dir, ss_lower_bound=-0.1, ss_upper_bound=1.0,
                    resolution=None):
    """ Draw a silhouette plot.
    Args:
        labels (array): array-like clustering results, {sample index: cluster label}
        frame_dr (ndarray): matrix after dimension reduction
        num_clusters (int): number of clusters in labels
        silhouette_avg (float): silhouette score
        out_dir (dir): path to output folder
        ss_lower_bound (float): The silhouette coefficient can range from -1, 1. This parameter sets the lower
        limit for plotting.
        ss_upper_bound (float): The silhouette coefficient can range from -1, 1. This parameter sets the upper
        limit for plotting.
        resolution (float): Louvain clustering resolution
    Returns:
        PDF image of silhouette plot
    """
    # Compute the silhouette scores for each cell
    sample_silhouette_values = silhouette_samples(frame_dr, labels)
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 7)
    y_lower = 10
    if min(labels) == 0:     # cluster labels must start from 1
        labels = labels + 1
    colors = plt.cm.jet(np.linspace(0, 1, num_clusters))
    for i in range(1, num_clusters + 1):
        # Aggregate the silhouette scores for cells belonging to cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        # The silhouette coefficient can range from -1, 1
        ax.set_xlim([ss_lower_bound, ss_upper_bound])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax.set_ylim([0, len(frame_dr) + (num_clusters + 1) * 10])

        # color = cm.nipy_spectral(float(i) / num_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                         0, ith_cluster_silhouette_values,
                         facecolor=colors[i-1], edgecolor=colors[i-1], alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks(np.arange(ss_lower_bound, ss_upper_bound+0.1, 0.2))
    if resolution:
        out_pdf_file = '{}/silhouette_{}_{}_{}.pdf'.format(out_dir, frame_dr.shape[1], num_clusters, resolution)
    else:
        out_pdf_file = '{}/silhouette_{}_{}.pdf'.format(out_dir, frame_dr.shape[1], num_clusters)
    plt.savefig(out_pdf_file, bbox_inches="tight")
    return sample_silhouette_values

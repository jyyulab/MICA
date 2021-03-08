#!/usr/bin/env python3

import sys
import numpy as np
import anndata
import logging
import networkx as nx
import time
import argparse
import fast_histogram
import community
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE
from functools import partial
import umap


def main():
    head_description = '''MICA is a mutual information-based nonlinear clustering analysis tool designed for 
    scRNA-seq data.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)
    parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                        help='Path to an input file (h5ad file or tab-delimited text file)')
    parser.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                        help='Path to final output directory')
    parser.add_argument('-d', '--dr-dim', metavar='INT', required=False, default=2, type=int,
                        help='Number of dimensions to reduce features of the input matrix to')
    parser.add_argument('-m', '--dist-metric', metavar='STR', required=False, default="mi", type=str,
                        help='Distance metric for UMAP dimension reduction [mi | euclidean] (default: mi)')
    parser.add_argument('-n', '--num-neighbors', metavar='INT', required=False, default=30, type=int,
                        help='Number of neighbors of building neighboring graph (default: 20)')
    parser.add_argument('-e', '--resolution', metavar='FLOAT', required=False, default=1.0, type=float,
                        help='Determines size of the communities. (default: 1.0)')
    parser.add_argument('-j', '--num-jobs', metavar='INT', required=False, default=1, type=int,
                        help='Number of jobs in building the neighbor graph in parallel (default: 1)')
    parser.add_argument('-v', '--visual-method', metavar='STR', required=False, default='umap', type=str,
                        help='Visualization embedding method [umap | tsne] (default: umap)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    frame = read_preprocessed_mat(args.input_file)
    # print(frame.iloc[500:4600,])
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Dimension reduction to {} dimensions using UMAP method and {} distance ...'.format(
        args.dr_dim, args.dist_metric))
    # frame = frame.iloc[500:4600,]
    frame_dr = dim_reduce(frame, dist=args.dist_metric, dim=args.dr_dim, out_dir=args.output_dir)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building {}-based {}-nearest neighbor graph ...'.format('euclidean', args.num_neighbors))
    G = build_graph(frame_dr, num_jobs=args.num_jobs)
    nx.write_gml(G, '{}/knn_{}_{}.gml'.format(args.output_dir, 'euclidean', args.num_neighbors))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Performing graph clustering ...')
    partition = graph_clustering(G)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Visualizing clustering results using {} ...'.format(args.visual_method))
    visual_embed(partition, frame.index, frame_dr, args.output_dir, embed_method=args.visual_method)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


def read_preprocessed_mat(in_file):
    """Read in preprocessed matrix file into a dataframe."""
    if in_file.endswith('.txt'):
        frame = pd.read_csv(in_file, sep="\t", index_col=0).iloc[:, 0:]
    if in_file.endswith('.h5ad') or in_file.endswith('.h5'):
        adata = anndata.read_h5ad(in_file)
        frame = adata.to_df()
    return frame


def calc_norm_mi(arr1, arr2, bins, m):
    """ Calculates a normalized mutual information distance D(X, Y) = 1 - I(X, Y)/H(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (pandas series): gene expression data for cell 1
        arr2 (pandas series): gene expression data for cell 2
        marginals  (ndarray): marginal probability matrix
        index1         (int): index of cell 1
        index2         (int): index of cell 2
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        a float between 0 and 1
    """
    fq = fast_histogram.histogram2d(arr1, arr2, range=[[arr1.min(), arr1.max()+1e-9], [arr2.min(), arr2.max()+1e-9]],
                                    bins=(bins, bins)) / float(m)
    sm = np.sum(fq * float(m), axis=1)
    tm = np.sum(fq * float(m), axis=0)
    sm = np.asmatrix(sm / float(sm.sum()))
    tm = np.asmatrix(tm / float(tm.sum()))
    sm_tm = np.matmul(np.transpose(sm), tm)
    div = np.divide(fq, sm_tm, where=sm_tm != 0, out=np.zeros_like(fq))
    ent = np.log(div, where=div != 0, out=np.zeros_like(div))
    agg = np.multiply(fq, ent, out=np.zeros_like(fq), where=fq != 0)
    joint_ent = -np.multiply(fq, np.log(fq, where=fq != 0, out=np.zeros_like(fq)),
                             out=np.zeros_like(fq), where=fq != 0).sum()
    return (joint_ent - agg.sum()) / joint_ent


def dim_reduce(frame, dim=2, n_neighbors=30, dist='mi', out_dir=None):
    """ Dimension reduction on a n_obs * n_vars matrix.
    Args:
        df (dataframe): preprocessed expression matrix as a dataframe
        dim (int): dimension to reduce n_vars to
        out_dir (str): path to output directory
    Returns:
        n_obs * dim numpy.ndarray
    """
    if dist == 'euclidean':
        clusterable_embedding = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.0,
                                          n_components=dim).fit_transform(frame.to_numpy())
    elif dist == 'mi':
        num_bins = int((frame.shape[0]) ** (1 / 3.0))
        num_genes = frame.shape[1]
        metric_params = {"bins": num_bins, "m": num_genes}
        clusterable_embedding = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.0,
                                          n_components=dim, metric=calc_norm_mi, metric_kwds=metric_params
                                          ).fit_transform(frame.to_numpy())
    else:
        sys.exit('Error - distance metric not supported: {}'.format(dist))
    if out_dir:
        np.savetxt('{}/mat_dr.csv'.format(out_dir), clusterable_embedding, delimiter=',')
    return clusterable_embedding


def build_graph(frame_dr, dis_metric='euclidean', num_neighbors=30, knn_algorithm='ball_tree', num_jobs=1):
    """Build a graph representation of the dimension reduced matrix.
    Args:
        frame_dr (numpy ndarray): dimension reduced n_obs * dim matrix
        dis_metric (str): 'MI' or 'euclidean'
        num_neighbors (int): number of neighbors
        knn_algorithm (str): algorithm used to compute the nearest neighbors
        num_jobs (None or int): n_jobs parameter in sklearn.neighbors.NearestNeighbors
    Returns:
        networkx graph
    """
    if dis_metric == 'mi':
        num_bins = int((frame_dr.shape[0]) ** (1 / 3.0))
        num_genes = frame_dr.shape[1]
        metric_params = {"bins": num_bins, "m": num_genes}
        nbrs = NearestNeighbors(n_neighbors=num_neighbors, algorithm=knn_algorithm, metric=calc_norm_mi,
                                metric_params=metric_params, n_jobs=num_jobs)
    elif dis_metric == 'euclidean':
        nbrs = NearestNeighbors(n_neighbors=num_neighbors, algorithm=knn_algorithm, n_jobs=num_jobs)
    else:
        sys.exit('Error - invalid distance metric: {}'.format(dis_metric))
    nbrs.fit(frame_dr)
    logging.info(nbrs.get_params())
    kneighbor_graph = nbrs.kneighbors_graph(frame_dr, mode='distance').toarray()
    return nx.from_numpy_matrix(kneighbor_graph)


def graph_clustering(G, method='louvein', resolution=1.0):
    """ Perform graph-based clustering.
    Args:
        G (nx graph): G to perform community detection
        method (str): louvein
        resolution (float): Determines size of the communities. (default: 1.0)
    Returns:
    """
    if method == 'louvein':
        partition = community.best_partition(G, resolution=resolution)
        logging.info('Clustering labels: {}'.format(set(partition.values())))
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    return partition


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


def visual_embed(partition, index, frame_dr, out_dir, embed_method='umap', dis_metric='euclidean', perplexity=30):
    """ Visualize clustering results using tSNE.
    Args:
        partition (dict): clustering results, {index: cluster label}
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
        partial_calc_norm_mi = partial(calc_norm_mi, bins=num_bins, m=num_genes)
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
    out_txt_file = '{}/clustering_{}_{}_{}.txt'.format(out_dir, embed_method, dis_metric, frame_dr.shape[1])
    final_res.to_csv(out_txt_file, sep="\t")
    out_png_file = '{}/clustering_{}_{}_{}.png'.format(out_dir, embed_method, dis_metric, frame_dr.shape[1])
    scatter_plot(final_res, out_png_file)


if __name__ == "__main__":
    main()

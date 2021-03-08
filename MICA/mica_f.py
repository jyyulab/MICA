#!/usr/bin/env python3

import sys
import numpy as np
import anndata
import logging
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
import networkx as nx
import time
import argparse
import fast_histogram
import community
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
    parser.add_argument('-r', '--dr-method', metavar='STR', required=False, default='pca', type=str,
                        help='Dimension reduction method [pca | mds | none] (default: pca)')
    parser.add_argument('-d', '--dr-dim', metavar='INT', required=False, default=100, type=int,
                        help='Number of dimensions to reduce features of the input matrix to')
    parser.add_argument('-m', '--dist-metric', metavar='STR', required=False, default="mi", type=str,
                        help='Distance metric for calculating neighboring cells [mi | euclidean] (default: mi)')
    parser.add_argument('-n', '--num-neighbors', metavar='INT', required=False, default=20, type=int,
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
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    if args.dr_method == 'none':
        logging.info('Skipping dimension reduction ... ')
        frame_dr = frame.to_numpy()
    else:
        start = time.time()
        logging.info('Dimension reduction to {} dimensions using {} method ...'.format(args.dr_dim, args.dr_method))
        frame_dr = dim_reduce(frame, dim=args.dr_dim, out_dir=args.output_dir)
        end = time.time()
        runtime = end - start
        logging.info('Done. Runtime: {} seconds'.format(runtime))

    # start = time.time()
    # logging.info('Calculating cell-by-cell distance matrix using {} method ...'.format(args.dist_metric))
    # df_dr = pd.DataFrame(frame_dr, index=frame.index)
    # num_bins = int((frame_dr.shape[0]) ** (1 / 3.0))
    # num_genes = frame_dr.shape[1]
    # dist_mat = calc_dis_mat(df_dr, df_dr, num_bins, num_genes)
    # dist_mat.to_csv('{}/dist_mat_{}.csv'.format(args.output_dir, args.dist_metric))
    # end = time.time()
    # runtime = end - start
    # logging.info('Done. Runtime: {} seconds'.format(runtime))

    # start = time.time()
    # logging.info('Calculating cell-by-cell distance matrix using {} method ...'.format(args.dist_metric))
    # num_bins = int((frame_dr.shape[0]) ** (1 / 3.0))
    # num_genes = frame_dr.shape[1]
    # dist_mat = calc_dis_mat_np(frame_dr, frame_dr, num_bins, num_genes, frame.index)
    # dist_mat.to_csv('{}/dist_mat_np_{}_{}.csv'.format(args.output_dir, args.dist_metric, args.dr_dim))
    # end = time.time()
    # runtime = end - start
    # logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building {}-based {}-nearest neighbor graph ...'.format(args.dist_metric, args.num_neighbors))
    G = build_graph(frame_dr, dis_metric=args.dist_metric, num_jobs=args.num_jobs)
    nx.write_gml(G, '{}/knn_{}_{}.gml'.format(args.output_dir, args.dist_metric, args.num_neighbors))
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
    visual_embed(partition, frame.index, frame_dr, args.output_dir, embed_method=args.visual_method,
                 dis_metric=args.dist_metric)
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


def calc_mi_f(arr1, arr2, bins, m):
    """ Calculates bin-based calculation for mutual information (MI) in between two arrays. It builds a 2d histogram,
    which is used to calculate P(arr1, arr2)
    Args:
        arr1 (pandas series): gene expression data for a given cell
        arr2 (pandas series): gene expression for another cell
        bins           (int): number of bins
        m              (int): number of genes
    Return:
        MI value
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
    return agg.sum()


def calc_dis_mat(mat1, mat2, bins, m):
    """ Wrapper of calc_mi for calculating mutual information for two matrices
    Args:
        mat1 (pandas dataframe): exp matrix of a slice of cells, with cells as rows from original file
                                  and all gene expression attributes as columns
        mat2 (pandas dataframe): exp matrix of another slice of cells
        bins              (int): number of bins
        m                 (int): number of genes
    Returns:
        df (pandas dataframe with dimension mat1.index * mat2.index)
    """
    df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
    for c in mat2.index:
        df.loc[mat1.index, c] = mat1.apply(calc_norm_mi, axis=1, args=(mat2.loc[c, :], bins, m))
    return df


def calc_dis_mat_np(mat1, mat2, bins, m, index):
    """ Wrapper of calc_mi for calculating mutual information for two matrices(numpy.ndarray)
    Args:
        mat1 (numpy.ndarray): exp matrix of a slice of cells, with cells as rows from original file
                                  and all gene expression attributes as columns
        mat2 (numpy.ndarray): exp matrix of another slice of cells
        bins           (int): number of bins
        m              (int): number of genes
        index       (series): cell index
    Returns:
        df (pandas dataframe with dimension mat1.index * mat2.index)
    """
    df = pd.DataFrame(data=0, index=index, columns=index)
    for i, c in enumerate(index):
        df.loc[index, c] = np.apply_along_axis(calc_norm_mi, 1, mat1, mat2[i, :], bins, m)
    return df


def dim_reduce(df, dim=100, method='pca', num_jobs=None, out_dir=None):
    """ Dimension reduction on a n_obs * n_vars matrix.
    Args:
        df (dataframe): preprocessed expression matrix as a dataframe
        dim (int): dimension to reduce n_vars to
        method (str): PCA or MDS
        num_jobs (None or int): n_jobs parameter in sklearn.manifold.MDS
        out_dir (str): path to output directory
    Returns:
        n_obs * dim numpy.ndarray
    """
    if method == 'pca':
        embedding = PCA(n_components=dim)
    elif method == 'mds':
        embedding = MDS(n_components=dim, n_jobs=num_jobs)
    else:
        sys.exit('Error - invalid dimension reduction method: {}'.format(method))
    frame_dr = embedding.fit_transform(df)
    if out_dir:
        np.savetxt('{}/mat_dr.csv'.format(out_dir), frame_dr, delimiter=',')
    return frame_dr


def build_graph(frame_dr, dis_metric='mi', num_neighbors=20, knn_algorithm='ball_tree', num_jobs=1):
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


def visual_embed(partition, index, frame_dr, out_dir, embed_method='umap', dis_metric='mi', perplexity=30):
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
            res = umap.UMAP(random_state=30, metric=partial_calc_norm_mi, n_neighbors=20,
                            min_dist=0.25).fit_transform(frame_dr)
            embed = pd.DataFrame(data=res, index=index, columns=["X", "Y"])
        else:
            sys.exit('Error - invalid embed method: {}'.format(embed_method))
    else:
        if embed_method == 'tsne':
            embed = TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=perplexity, random_state=10,
                         early_exaggeration=12.0).fit_transform(frame_dr)
        elif embed_method == 'umap':
            res = umap.UMAP(random_state=30, metric='euclidean', n_neighbors=20, min_dist=0.25).fit_transform(frame_dr)
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

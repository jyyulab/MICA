#!/usr/bin/env python3

import sys
import numpy as np
import umap
import logging
import node2vec
from .distance import numba_calc_mi_dis
from .aux_utils import run_shell_command


def dim_reduce_umap(frame, dim=12, n_neighbors=30, dist='mi', out_dir=None):
    """ Dimension reduction on a n_obs * n_vars matrix.
    Args:
        frame (data frame): preprocessed expression matrix as a dataframe
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
                                          n_components=dim, metric=numba_calc_mi_dis, metric_kwds=metric_params
                                          ).fit_transform(frame.to_numpy())
    else:
        sys.exit('Error - distance metric not supported: {}'.format(dist))
    if out_dir:
        np.savetxt('{}/mat_dr.csv'.format(out_dir), clusterable_embedding, delimiter=',')
    return clusterable_embedding


def dim_reduce_node2vec(graph, dim=12, walk_len=30, n_walks=40, n_workers=1):
    """ Dimension reduction on a n_obs * n_vars matrix using node2vec python implementation.
    Args:
        graph (Networkx Graph)
        dim (int): dimension to reduce nodes to (default: 12)
        walk_len (str): path to output directory (default: 30)
        n_walks (str): number of random walks per node (default: 40)
        n_workers (str): Number of workers for parallel execution (default: 1)
    Returns:

    """
    n2v = node2vec.Node2Vec(graph, dimensions=dim, walk_length=walk_len, num_walks=n_walks, workers=n_workers,
                            weight_key='MI')
    model = n2v.fit(window=10, min_count=1, batch_words=4)
    return model.wv


def dim_reduce_node2vec_hp(edgelist_file, out_emb_file, dim=12, walk_len=60, n_walks=120):
    """ Dimension reduction on a n_obs * n_vars matrix using node2vec C++ implementation.
    Args:
        edgelist_file (txt file): path to a graph edgelist file,
                                  (format: node1_id_int node2_id_int <weight_float, optional>)
        out_emb_file (txt file): path ot output embedding file
        dim (int): dimension to reduce nodes to (default: 12)
        walk_len (str): path to output directory (default: 30)
        n_walks (str): number of random walks per node (default: 40)
    Returns:
        out_emb_file (txt file): output embedding file
    """
    cmd = 'node2vec -i:{} -o:{} -d:{} -l:{} -r:{} -w'.format(edgelist_file, out_emb_file, dim, walk_len, n_walks)
    logging.info(cmd)
    run_shell_command(cmd)


def dim_reduce_deepwalk(edgelist_file, out_emb_file, dim=12, walk_len=60, n_walks=120):
    """ Dimension reduction on a n_obs * n_vars matrix.
    Args:
        edgelist_file (txt file): path to a graph edgelist file,
                                  (format: node1_id_int node2_id_int <weight_float, optional>)
        out_emb_file (txt file): path ot output embedding file
        dim (int): dimension to reduce nodes to (default: 10)
        walk_len (str): path to output directory (default: 30)
        n_walks (str): number of random walks per node (default: 200)
    Returns:
        out_emb_file (txt file): output embedding file
    """
    cmd = 'deepwalk --input {} --format {} --output {} --representation-size {} --number-walks {} ' \
          '--walk-length {}'.format(edgelist_file, format, out_emb_file, dim, n_walks, walk_len)
    logging.info(cmd)
    run_shell_command(cmd)

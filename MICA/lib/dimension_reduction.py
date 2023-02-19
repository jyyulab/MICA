#!/usr/bin/env python3

import sys
import numpy as np
import umap
import logging
# import node2vec
import torch
import torch_geometric
from torch_geometric.nn.models import Node2Vec
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from .distance import numba_calc_mi_dis
from .aux_utils import run_shell_command


def dim_reduce_global(df, dim=50, method='pca', num_jobs=None, out_dir=None):
    """ Dimension reduction on a n_obs * n_vars matrix.
    Args:
        df (ndarray): preprocessed expression matrix as a dataframe
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
        num_bins = int((frame.shape[1]) ** (1 / 3.0))
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


# def dim_reduce_node2vec(graph, out_emb_file, dim=12, walk_len=60, n_walks=120, n_workers=1):
    """ Dimension reduction on a n_obs * n_vars matrix using node2vec python implementation.
    Args:
        graph (Networkx Graph)
        out_emb_file (txt file): path ot output embedding file
        dim (int): dimension to reduce nodes to (default: 12)
        walk_len (str): path to output directory (default: 60)
        n_walks (str): number of random walks per node (default: 120)
        n_workers (str): Number of workers for parallel execution (default: 1)
    Returns:
        out_emb_file (txt file): output embedding file
    """
    # n2v = node2vec.Node2Vec(graph, dimensions=dim, walk_length=walk_len, num_walks=n_walks, workers=n_workers,
    #                         weight_key='MI')
    # model = n2v.fit(window=10, min_count=1, batch_words=4)
    # model.wv.save_word2vec_format(out_emb_file)


def dim_reduce_node2vec_pecanpy(edgelist_file, out_emb_file, mode='SparseOTF', dim=20, walk_len=100, n_walks=120,
                                context_size=20, num_jobs=10, hyper_p=0.5, hyper_q=0.5):
    """ Dimension reduction on a n_obs * n_vars matrix using node2vec pecanpy implementation.
    Args:
        edgelist_file (txt file): path to a graph edgelist file,
                                  (format: node1_id_int node2_id_int <weight_float, optional>)
        mode (str): PreComp or SparseOTF or DenseOTF (default: SparseOTF)
        out_emb_file (txt file): path ot output embedding file
        dim (int): dimension to reduce nodes to (default: 20)
        walk_len (str): Length of walk per source (default: 100)
        n_walks (str): number of random walks per node (default: 120)
        context_size (int): window size of a random walk. (default: 20)
        num_jobs (int): Number of parallel workers (default: 10)
        hyper_p (float): Return hyperparameter controls the probability of a walk staying inward revisiting nodes.
        (default: 0.5)
        hyper_q (float): Inout hyperparameter controls the probability of a walk staying close to the preceeding nodes
        or moving outward farther away. (default: 0.5)
    Returns:
        out_emb_file (txt file): output embedding file
    """
    cmd = 'pecanpy --input {} --output {} --mode {} --dimensions {} --walk-length {} ' \
          '--num-walks {} --window-size {} --workers {} --p {} --q {} --weighted'.format(edgelist_file, out_emb_file,
                                                                                         mode, dim, walk_len, n_walks,
                                                                                         context_size, num_jobs,
                                                                                         hyper_p, hyper_q)
    logging.info(cmd)
    run_shell_command(cmd)


def dim_reduce_node2vec_pytorch_geometric(knn_graph, out_emb_file, sparse=True, dim=20, walk_len=100, n_walks=120,
                                          context_size=20, hyper_p=0.5, hyper_q=0.5):
    graph_data = torch_geometric.utils.from_networkx(knn_graph)
    pyg_model = Node2Vec(graph_data.edge_index, embedding_dim=dim, walk_length=walk_len,
                         context_size=context_size, walks_per_node=n_walks, p=hyper_p,
                         q=hyper_q, sparse=sparse).to("cuda" if torch.cuda.is_available() else "cpu")
    pyg_forward = pyg_model.forward().detach().numpy()

    with open(out_emb_file, "w") as f:
        f.write(str(pyg_forward.shape[0]) + " " + str(pyg_forward.shape[1]))
        f.write("\n")
        for i in range(len(pyg_forward)):
            w = pyg_forward[i]
            f.write(str(i))
            f.write(" ")
            f.write(" ".join(str(e) for e in w))
            f.write("\n")


def dim_reduce_node2vec_c(edgelist_file, out_emb_file, dim=12, walk_len=60, n_walks=120):
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
    format = 'edgelist'
    cmd = 'deepwalk --input {} --format {} --output {} --representation-size {} --number-walks {} ' \
          '--walk-length {}'.format(edgelist_file, format, out_emb_file, dim, n_walks, walk_len)
    logging.info(cmd)
    run_shell_command(cmd)

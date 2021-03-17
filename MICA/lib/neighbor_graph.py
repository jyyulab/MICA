#!/usr/bin/env python3

import sys
import logging
import networkx as nx
from sklearn.neighbors import NearestNeighbors
import umap
from .distance import numba_calc_mi_dis
from .distance import calc_norm_mi


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
        nbrs = NearestNeighbors(n_neighbors=num_neighbors, algorithm=knn_algorithm, metric=numba_calc_mi_dis,
                                metric_params=metric_params, n_jobs=num_jobs)
    elif dis_metric == 'euclidean':
        nbrs = NearestNeighbors(n_neighbors=num_neighbors, algorithm=knn_algorithm, n_jobs=num_jobs)
    else:
        sys.exit('Error - invalid distance metric: {}'.format(dis_metric))
    nbrs.fit(frame_dr)
    logging.info(nbrs.get_params())
    kneighbor_graph = nbrs.kneighbors_graph(frame_dr, mode='distance').toarray()
    return nx.from_numpy_matrix(kneighbor_graph)


def nearest_neighbors_umap(mat, num_neighbors=30, angular=False, random_state=None):
    """ Build a graph representation of the dimension reduced matrix.
    Args:
        mat (numpy ndarray): n_obs * n_var matrix
        num_neighbors (int): number of neighbors
        angular (bool): whether to use angular rp trees in NN approximation (default: False)
        random_state (np.random): the random state to use for approximate NN computations (default: None)
            If int, random_state is the seed used by the random number generator;
            If RandomState instance, random_state is the random number generator;
            If None, the random number generator is the RandomState instance used
            by `np.random`.
    Returns:
        knn_indices: array of shape (n_samples, n_neighbors)
            The indices on the ``n_neighbors`` closest points in the dataset.
        knn_dists: array of shape (n_samples, n_neighbors)
            The distances to the ``n_neighbors`` closest points in the dataset.
        rp_forest: list of trees
            The random projection forest used for searching (if used, None otherwise)
    """
    num_bins = int((mat.shape[0]) ** (1 / 3.0))
    num_genes = mat.shape[1]
    metric_params = {"bins": num_bins, "m": num_genes}
    knn_indices, knn_dists, forest = umap.umap_.nearest_neighbors(X=mat, n_neighbors=num_neighbors,
                                                                  metric=numba_calc_mi_dis, metric_kwds=metric_params,
                                                                  angular=angular, random_state=random_state)
    return knn_indices, knn_dists, forest


def build_graph_from_indices(knn_indices, knn_dists):
    """ Build a NetworkX undirected graph from nearest neighbor indices and distances.
    Args:
        knn_indices: array of shape (n_samples, n_neighbors)
            The indices on the ``n_neighbors`` closest points in the dataset.
        knn_dists: array of shape (n_samples, n_neighbors)
            The distances to the ``n_neighbors`` closest points in the dataset.
    Returns:
        A Networkx undirected graph
    """
    knn_graph = nx.Graph()
    for c_index, cell in enumerate(knn_indices):
        knn_graph.add_node(cell[0])
        for n_index, neighbor in enumerate(cell):
            if n_index == 0:
                continue
            knn_graph.add_node(neighbor)
            knn_graph.add_edge(cell[0], neighbor, MI=knn_dists[c_index][n_index])
    return knn_graph

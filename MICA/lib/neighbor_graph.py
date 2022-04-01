#!/usr/bin/env python3

import sys
import logging
import networkx as nx
from sklearn.neighbors import NearestNeighbors
from pynndescent import NNDescent
from .distance import numba_calc_mi_dis


def build_graph(frame_dr, dis_metric='euclidean', num_neighbors=20, knn_algorithm='ball_tree', num_jobs=1):
    """ Build a graph representation of the dimension reduced matrix.
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
        # Use n^{1/3} as the bin size, where n is the number of genes.
        num_bins = int((frame_dr.shape[1]) ** (1 / 3.0))
        logging.info('Number of bins for estimating MI: {}'.format(num_bins))
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


def nearest_neighbors_NNDescent(mat, num_neighbors=100, pruning_degree_multi=3.0, diversify_p=0.0, num_jobs=10):
    """ Build a graph representation of the dimension reduced matrix.
    Args:
        mat (numpy ndarray): n_obs * n_var matrix
        num_neighbors: int (optional, default=30)
            The number of neighbors to use in k-neighbor graph graph_data structure
            used for fast approximate nearest neighbor search. Larger values
            will result in more accurate search results at the cost of
            computation time.
        pruning_degree_multi: float (optional, default=3.0)
            How aggressively to prune the graph. Since the search graph is undirected
            (and thus includes nearest neighbors and reverse nearest neighbors) vertices
            can have very high degree -- the graph will be pruned such that no
            vertex has degree greater than pruning_degree_multiplier * n_neighbors.
        diversify_p: float (optional, default=1.0)
            The search graph get "diversified" by removing potentially unnecessary
            edges. This controls the volume of edges removed. A value of 0.0 ensures
            that no edges get removed, and larger values result in significantly more
            aggressive edge removal. A value of 1.0 will prune all edges that it can.
        num_jobs: int or None, optional (default=None)
            The number of parallel jobs to run for neighbors index construction.
            ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
            ``-1`` means using all processors.
    Returns:
        knn_indices: array of shape (n_samples, n_neighbors)
            The indices on the ``n_neighbors`` closest points in the dataset.
        knn_dists: array of shape (n_samples, n_neighbors)
            The distances to the ``n_neighbors`` closest points in the dataset.
    """
    # Use n^{1/3} as the bin size, where n is the number of genes.
    num_bins = int((mat.shape[1]) ** (1 / 3.0))
    logging.info('Number of bins for estimating MI: {}'.format(num_bins))
    num_genes = mat.shape[1]
    metric_params = {"bins": num_bins, "m": num_genes}
    knn_indices, knn_dists = NNDescent(mat, n_neighbors=num_neighbors, metric=numba_calc_mi_dis,
                                       metric_kwds=metric_params, pruning_degree_multiplier=pruning_degree_multi,
                                       diversify_prob=diversify_p, n_jobs=num_jobs).neighbor_graph
    return knn_indices, knn_dists


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

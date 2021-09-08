#!/usr/bin/env python3

import sys
import numpy as np
import community
import logging
from multiprocessing import Pool
from functools import partial
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score
from MICA.lib import visualize as vs


def graph_clustering(G, method='louvain', min_resolution=0.2, max_resolution=3.4, step_size=0.4):
    """ Perform graph-based clustering.
    Args:
        G (nx graph): G to perform community detection
        method (str): Louvain algorithm
        min_resolution (float): Determines minimum size of the communities. (default: 0.2)
        max_resolution (float): Determines maximum size of the communities. (default: 3.4)

    Returns:
        Clustering results
    """
    partitions = []
    if method == 'louvain':
        for resolution in np.arange(min_resolution, max_resolution+0.1, step_size):
            partition = community.best_partition(G, resolution=resolution)
            # logging.info('Clustering labels: {}'.format(set(partition.values())))
            partitions.append(partition)
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    return partitions


def best_partition_wrapper(G, max_resolution):
    return community.best_partition(G, resolution=max_resolution)


def graph_clustering_parallel(G, method='louvain', min_resolution=0.2, max_resolution=3.4, step_size=0.4,
                              num_workers=10):
    """ Perform graph-based clustering in parallel.
    Args:
        G (nx graph): G to perform community detection
        method (str): Louvain algorithm
        min_resolution (float): Determines minimum size of the communities. (default: 0.4)
        max_resolution (float): Determines maximum size of the communities. (default: 3.4)
        step_size (float): step size to sweep resolution from min_resolution to max_resolution
        num_workers (int): number of processes (default: 10)
    Returns:
        Clustering results
    """
    pool = Pool(processes=num_workers)
    if method == 'louvain':
        community_partial = partial(best_partition_wrapper, G)
        partitions = pool.map(community_partial, list(np.arange(min_resolution, max_resolution+0.1, step_size)))
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    pool.close()
    return partitions


def silhouette_analysis(clustering_res, num_clusters, frame_dr, out_dir):
    """ Calculate silhouette scores and draw silhouette plot of consensus clustering results for a given number of
    clusters.
    Args:
        clustering_res (dict): clustering results, {sample index: cluster label}
        num_clusters (int): number of clusters
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
    Outputs:
        silhouette score (float)
        PDF image of silhouette plot
    """
    if num_clusters == 1:
        # Number of clusters is 1. Set silhouette score to 0.0
        logging.info('Number of clusters: {}, silhouette_avg: {}'.format(num_clusters, 0.0))
        return 0.0
    labels = clustering_res['label']
    silhouette_avg = silhouette_score(frame_dr, labels)
    logging.info('Number of clusters: {}, silhouette score: {}'.format(num_clusters, silhouette_avg))
    vs.silhouette_plot(labels, frame_dr, num_clusters, silhouette_avg, out_dir)
    return silhouette_avg

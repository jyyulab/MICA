#!/usr/bin/env python3

import sys
import math
import numpy as np
import community
import logging
from multiprocessing import Pool
from functools import partial
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'Arial'
from sklearn.metrics import silhouette_score
from MICA.lib import visualize as vs


def graph_clustering(G, method='louvain', min_resolution=-2.0, max_resolution=3.0, step_size=0.2):
    """ Perform graph-based clustering.
    Args:
        G (nx graph): G to perform community detection
        method (str): Louvain algorithm
        min_resolution (float): Determines minimum size of the communities. (default: -2.0)
        max_resolution (float): Determines maximum size of the communities. (default: 3.0)
        step_size (float): Determines the step size for sweeping the resolutions
    Returns:
        Clustering results
    """
    partitions = []
    if method == 'louvain':
        reso_lst = []
        for reso in list(np.arange(min_resolution, max_resolution + 0.1, step_size)):
            reso_lst.append(math.exp(reso))
        for resolution in reso_lst:
            partition = community.best_partition(G, resolution=resolution)
            # logging.info('Clustering labels: {}'.format(set(partition.values())))
            partitions.append(partition)
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    return partitions


def best_partition_wrapper(G, max_resolution):
    return community.best_partition(G, resolution=max_resolution), max_resolution


def graph_clustering_parallel(G, method='louvain', min_resolution=-2.0, max_resolution=3.0, step_size=0.2,
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
        partition_resolutions [(dict, float)]: a list of clustering_result and resolution tuples
    """
    pool = Pool(processes=num_workers)
    reso_lst = []
    for reso in list(np.arange(min_resolution, max_resolution+0.001, step_size)):
        reso_lst.append(math.exp(reso))
    if method == 'louvain':
        community_partial = partial(best_partition_wrapper, G)
        partition_resolutions = pool.map(community_partial, reso_lst)
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    pool.close()
    return partition_resolutions


def silhouette_analysis(clustering_res, num_clusters, frame_dr, out_dir, resolution=None):
    """ Calculate silhouette scores and draw silhouette plot of consensus clustering results for a given number of
    clusters.
    Args:
        clustering_res (dict): clustering results, {sample index: cluster label}
        num_clusters (int): number of clusters
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
        resolution (float): Louvain clustering resolution
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
    vs.silhouette_plot(labels, frame_dr, num_clusters, silhouette_avg, out_dir, ss_lower_bound=-0.6, ss_upper_bound=1.0,
                       resolution=resolution)
    return silhouette_avg

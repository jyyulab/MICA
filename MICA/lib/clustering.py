#!/usr/bin/env python3

import sys
import numpy as np
import community
from multiprocessing import Pool
from functools import partial


def graph_clustering(G, method='louvain', max_resolution=3.4):
    """ Perform graph-based clustering.
    Args:
        G (nx graph): G to perform community detection
        method (str): Louvain algorithm
        max_resolution (float): Determines maximum size of the communities. (default: 3.4)
    Returns:
        Clustering results
    """
    partitions = []
    if method == 'louvain':
        for resolution in range(1.0, max_resolution+0.1, 0.4):
            partition = community.best_partition(G, resolution=resolution)
            # logging.info('Clustering labels: {}'.format(set(partition.values())))
            partitions.append(partition)
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    return partitions


def best_partition_wrapper(G, max_resolution):
    return community.best_partition(G, resolution=max_resolution)


def graph_clustering_parallel(G, method='louvain', min_resolution=0.4, max_resolution=3.4, step_size=0.4,
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

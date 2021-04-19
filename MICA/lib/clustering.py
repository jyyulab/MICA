#!/usr/bin/env python3

import sys
import logging
import community


def graph_clustering(G, method='louvain', resolution=1.0):
    """ Perform graph-based clustering.
    Args:
        G (nx graph): G to perform community detection
        method (str): louvein
        resolution (float): Determines size of the communities. (default: 1.0)
    Returns:
    """
    if method == 'louvain':
        partition = community.best_partition(G, resolution=resolution)
        logging.info('Clustering labels: {}'.format(set(partition.values())))
    else:
        sys.exit('Error - invalid graph clustering method: {}'.format(method))
    return partition

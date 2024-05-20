#!/usr/bin/env python3

import numpy as np
import pandas as pd
from sklearn import cluster
import logging


def group_partition(partitions, index):
    """ Group partitions based on number of clusters.
    Args:
        partitions: each dataframe is a clustering results with two columns (cell_index, clustering_label)
    Returns:
        Clustering results
    """
    cluster_dict = dict()
    for partition in partitions:
        labels = [x + 1 for x in partition.values()]
        num_cluster = len(set(labels))
        clustering_res = pd.DataFrame(data=labels, index=index, columns=["label"])
        if num_cluster in cluster_dict.keys():
            cluster_dict[num_cluster].append(clustering_res)
        else:
            cluster_dict[num_cluster] = [clustering_res]
    return cluster_dict


def consensus_sc3(clustering_results, n_clusters, common_name=None):
    """ Implement SC3's consensus clustering. (https://www.nature.com/articles/nmeth.4236)
    Args:
        clustering_results (list of dataframes): each dataframe is a clustering results with two columns (cell_index,
                                                 clustering_label)
        n_clusters (int): number of clusters
        common_name (name): common string to name the output files
    Returns:
        Clustering results
    """
    n_iter = len(clustering_results)
    logging.info('Number of clusters: {}'.format(n_clusters))
    logging.info('Number of clustering results: {}'.format(n_iter))
    if n_iter == 0:
        return None
    if len(clustering_results) == 0:
        return None

    conss_binary_mat = np.zeros((clustering_results[0].shape[0], clustering_results[0].shape[0]))
    for i in range(n_iter):
        arr = clustering_results[i].to_numpy()
        mask = arr[:, None] == arr
        binary_mat = mask[:,:,0].astype(int)
        conss_binary_mat += binary_mat
    conss_binary_mat = conss_binary_mat / n_iter

    clust = cluster.AgglomerativeClustering(
        linkage="complete", n_clusters=n_clusters, affinity="euclidean"
    )
    clust.fit(conss_binary_mat)

    cclust = pd.DataFrame(data=clust.labels_, index=clustering_results[0].index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})

    out_file = common_name + "_k" + str(n_clusters)
    # out_file_hdf = out_file + "_cclust.h5"
    # cclust.to_hdf(out_file_hdf, "cclust")
    # conss_binary_mat.to_hdf(out_file_hdf, "membership")
    return cclust, out_file


def consensus_mlca(clustering_results, n_clusters):
    """ The meta-cLustering algorithm for clustering clusters. """
    n_iter = len(clustering_results)
    logging.info('Number of clusters: {}'.format(n_clusters))
    logging.info('Number of clustering results: {}'.format(n_iter))
    if n_iter == 0:
        return None
    if len(clustering_results) == 0:
        return None
    return


def purify():
    """ To do: purify mis-classified cells. """
    return

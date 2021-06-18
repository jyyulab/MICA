#!/usr/bin/env python3

import numpy as np
import pandas as pd
import logging
from sklearn import cluster


def consensus_sc3(km_results, n_clusters, common_name=None):
    """ Implement SC3's consensus clustering. (https://www.nature.com/articles/nmeth.4236)
    Args:
        km_results (list of dataframes): each dataframe is a clustering results with two columns (cell_index,
                                         clustering_label)
        n_clusters (int): number of clusters
        common_name (name): common string to name the output files
    Returns:
        Clustering results
    """
    n_iter = len(km_results)
    print('Number of k-mean results: {}'.format(n_iter))
    if n_iter == 0:
        return None
    if len(km_results) == 0:
        return None

    conss_binary_mat = np.zeros((km_results[0].shape[0], km_results[0].shape[0]))
    for i in range(n_iter):
        arr = km_results[i].to_numpy()
        mask = arr[:, None] == arr
        binary_mat = mask[:,:,0].astype(int)
        conss_binary_mat += binary_mat
    conss_binary_mat = conss_binary_mat / n_iter

    clust = cluster.AgglomerativeClustering(
        linkage="complete", n_clusters=n_clusters, affinity="euclidean"
    )
    clust.fit(conss_binary_mat)

    cclust = pd.DataFrame(data=clust.labels_, index=km_results[0].index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})

    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    # conss_binary_mat.to_hdf(out_file_hdf, "membership")
    print('consensus_sc3 is done')
    return cclust, out_file

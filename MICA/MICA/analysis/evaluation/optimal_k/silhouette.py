#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_samples, silhouette_score


def silhouette_analysis(clustering_res, num_cluster, frame_dr, out_dir):
    """ Calculate silhouette scores and draw silhouette plot of consensus clustering results for a given number of
    clusters.
    Args:
        clustering_res (dict): clustering results, {sample index: cluster label}
        num_cluster (int): number of clusters
        frame_dr (ndarray): matrix after dimension reduction
        out_dir (dir): path to output folder
    Outputs:
        silhouette score (float)
        PDF image of silhouette plot
    """
    labels = clustering_res['label']
    silhouette_avg = silhouette_score(frame_dr, labels)
    print('silhouette_avg: {}'.format(silhouette_avg))
    return silhouette_avg

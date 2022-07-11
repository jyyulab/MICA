#!/usr/bin/env python3

import numpy as np
import pandas as pd
import random


def create_metacells(clustering_results, G, cell_barcodes, out_dir, num_metacells=0):
    """ Create MetaCell membership file
    Args:
        clustering_results (dict): a number of clustering results with keys cells and values labels
        G (networkx graph): a kNN graph to select the nearest neighbors
        cell_barcodes (series): cell barcodes for cells in the clustering results
        out_dir (dir): path to output folder
        num_metacells (int): count of MetaCells to create
    Return:
        A MetaCell membership file
    """
    if num_metacells == 0:  # Use the clustering result with the minimum resolution
        cls_res = min(clustering_results, key=lambda x: x[1])[0]
        cls_res_df = pd.DataFrame(cls_res.items(), columns=['cell', 'label'], index=cell_barcodes)
        pick_cells(cls_res_df, G, out_dir)
    else:
        pass


def pick_cells(cls_res_df, G, out_dir):
    """ Pick cells to form metacells using the global kNN graph
    Args:
        cls_res_df (dataframe): a clustering result with two columns cell and label
        G (networkx graph): a kNN graph to select the nearest neighbors
        out_dir (dir): path to output folder
    Return:
        A MetaCell membership file
    """
    labels = set([x for x in cls_res_df['label']])
    with open('{}/metacell_membership.txt'.format(out_dir), 'w') as fout:
        fout.write('barcode\tcell_index\tclustering_label\n')
    for label in labels:
        cells = cls_res_df.loc[cls_res_df['label'] == label]
        pivot = random.choice(list(cells['cell']))
        # print('pivot: {}'.format(pivot))
        neighbors = np.array([n for n in G.neighbors(pivot)])
        num_neighbors = len(neighbors)
        weights = np.array([G[pivot][n]['weight'] for n in neighbors])
        idx = np.argpartition(weights, num_neighbors-1)
        # print(idx[:num_neighbors])
        # print(neighbors[idx[:num_neighbors]])
        # Overlap neighbors of the pivot cell with cells in the cluster
        selected_cells = set(neighbors[idx[:num_neighbors]]).intersection(set(cells['cell']))
        selected_cells.add(pivot)
        metacell_df = cells.loc[cells['cell'].isin(selected_cells)]
        metacell_df.to_csv('{}/metacell_membership.txt'.format(out_dir), mode='a', sep='\t', header=False)


def create_norm_sum_mat(raw_count_mat_file, mc_member_file):
    """ Create MetaCell matrix by normalizing cumulative summation of raw counts, normalize by the number of cells
        selected to create the MetaCell.
    Args:
        raw_count_mat_file (str): path to a raw count matrix file
        mc_member_file (str): path to a MetaCell membership file
    Return:
        a gene * metacell matrix file
    """
    # better to implement this function in R as MICA may not take a raw count matrix as input
    pass

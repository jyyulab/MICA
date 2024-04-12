#!/usr/bin/env python3
from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
import numpy as np
import os


def PBMC20k(dim, reso):
    true_label_file = '{}/datasets/SilverStd/PBMC_20k/PBMC_20K_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['type'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['type']]
    true_label['label'] = labels
    PBMC20k = '{}/datasets/SilverStd/PBMC_20k/PBMC_20K_preprocessed.h5ad'.format(mica_data_path)
    PBMC20k_out = '{}/outputs/SilverStd/PBMC_20K/{}_{}'.format(mica_data_path, dim, reso)
    # PBMC20k_out = '{}/outputs/PBMC_20k_no_ge/{}'.format(mica_data_path, reso)
    if not os.path.isdir(PBMC20k_out):
        os.makedirs(PBMC20k_out)
    cmd = 'mica -i {} -o {} -d {} -e {}'.format(PBMC20k, PBMC20k_out, dim, reso)
    # cmd = 'mica -i {} -o {} -e {}'.format(PBMC20k, PBMC20k_out, reso)
    print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    # predict_label_file = '{}/clustering_umap_euclidean_{}_{}.txt'.format(PBMC20k_out, dim, reso)
    # predict_label_file = '{}/clustering_umap_euclidean_21952_{}.txt'.format(PBMC20k_out, reso)
    # predict_label = pd.read_csv(predict_label_file, delimiter='\t')
    # return adjusted_rand_score(true_label['label'], predict_label['label'])


def create_cmds():
    for dim in range(8, 100, 4):
        for reso in np.arange(0.2, 10.0, 0.4):
            reso_round = np.round(reso, 1)
            PBMC20k(dim, reso_round)


if __name__ == "__main__":
    mica_data_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    # create_cmds()
    # summary()

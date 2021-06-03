#!/usr/bin/env python3
from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
import numpy as np
import os


def create_cmds():
    for dim in range(8, 100, 4):
        input_file = '{}/datasets/GoldernStd/kolod/Kolod_MICA_input.txt'.format(mica_data_path)
        out_dir = '{}/outputs/GoldernStd/Kolod/MICA_GE/{}'.format(mica_data_path, dim)

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        cmd = 'mica -i {} -o {} -d {} -s 0.0 -ir 0.4 -ar 2.0 -ss 0.2'.format(input_file, out_dir, dim)
        print(cmd)


def summary():
    summary_file = '{}/outputs/summary.txt'.format(mica_data_path)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/GoldernStd/kolod/kolod_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    # print(true_label)

    with open('{}/outputs/summary.txt'.format(mica_data_path), 'w') as fout:
        for dim in range(8, 100, 4):
            yan_out = '{}/outputs/GoldernStd/Kolod/MICA_GE/{}'.format(mica_data_path, dim)
            print(dim)
            for reso in np.arange(0.4, 2.1, 0.2):
                reso_round = np.round(reso, 2)
                # print(reso_round)
                predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(yan_out, dim, reso_round)
                # print(predict_label_file)
                predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
                # new_index = [int(s.replace('V', '')) for s in predict_label.index]
                # predict_label.index = new_index
                # print(predict_label)
                merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                # print(merged)
                ari = adjusted_rand_score(merged['label_x'], merged['label_y'])

                fout.write('Kolod\t{}\t{}\t{}\n'.format(dim, reso_round, ari))
            # break


if __name__ == "__main__":
    mica_data_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    # create_cmds()
    summary()

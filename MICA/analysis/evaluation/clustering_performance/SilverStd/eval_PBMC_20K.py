#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
from sklearn.metrics.cluster import adjusted_rand_score


def create_cmds_ge():
    cmd_file = '{}/cmd_eval.sh'.format(test_path)
    with open(cmd_file, 'w') as fout:
        for dim in range(8, 100, 4):
            out_dir = '{}/dim_{}'.format(output_path, dim)
            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)
            cmd = 'mica -i {} -o {} -d {} -ar 10.0'.format(input_file, out_dir, dim)
            fout.write(cmd)
            fout.write('\n')
    print('Done')


def create_cmds_mds():
    pass


def create_cmds_auto():
    pass


def calc_ARIs_ge():
    summary_file = '{}/summary_ge_sqrt2.txt'.format(output_path)
    print(summary_file)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    # print(true_label)
    with open(summary_file, 'w') as fout:
        for dim in range(8, 100, 4):
            print(dim)
            clustering_out_dir = '{}/dim_{}'.format(output_path, dim)
            for reso in np.arange(0.4, 10.1, 0.4):
                reso_round = np.round(reso, 2)
                # print(reso_round)
                predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(clustering_out_dir, dim,
                                                                                     reso_round)
                # print(predict_label_file)
                predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
                # new_index = [int(s.replace('V', '')) for s in predict_label.index]
                # predict_label.index = new_index
                # print(predict_label)
                merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                # print(merged)
                ari = adjusted_rand_score(merged['type'], merged['label'])
                fout.write('{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ari))
                # break
            # break
    print('Done')


def calc_ARIs_mds(author):
    pass


def calc_ARIs_auto(author):
    pass


if __name__ == "__main__":
    root_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    level = 'SilverStd'
    author = 'PBMC_20k'
    data_path = '{}/datasets/{}/{}'.format(root_path, level, author)
    output_path = '{}/outputs/{}/{}'.format(root_path, level, author)
    test_path = '{}/tests/{}/{}'.format(root_path, level, author)
    input_file = '{}/{}_preprocessed.h5ad'.format(data_path, author)
    true_label_file = '{}/{}_true_label.txt'.format(data_path, author)
    # create_cmds_ge()
    calc_ARIs_ge()

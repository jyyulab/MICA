#!/usr/bin/env python3
from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
import numpy as np
import os


def create_cmds():
    for dim in range(8, 100, 4):
        input_file = '{}/datasets/GoldernStd/buettner/Buttner_MICA_input.txt'.format(mica_data_path)
        out_dir = '{}/outputs/GoldernStd/Buttner/MICA_GE/{}'.format(mica_data_path, dim)

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        cmd = 'mica -i {} -o {} -d {} -s 0.0 -ir 0.4 -ar 2.0 -ss 0.2'.format(input_file, out_dir, dim)
        print(cmd)


def summary():
    summary_file = '{}/outputs/summary.txt'.format(mica_data_path)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/GoldernStd/buettner/buettner_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    # print(true_label)

    with open('{}/outputs/summary.txt'.format(mica_data_path), 'w') as fout:
        for dim in range(8, 100, 4):
            yan_out = '{}/outputs/GoldernStd/Buttner/MICA_GE/{}'.format(mica_data_path, dim)
            # print(dim)
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

                fout.write('Buttner\t{}\t{}\t{}\n'.format(dim, reso_round, ari))
            # break


def calcu_ARI(author):
    root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'

    true_label_file = '{}/datasets/GoldernStd/{}/{}_true_label.txt'.format(root_dir, author, author)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    # print(true_label)

    out_dir = '{}/outputs/GoldernStd/{}/MICA_MDS/mica-e3e950ad-5991-4c68-adaf-3061519b1b82'.format(root_dir, author)
    cluster_mem_file = '{}/cwl_lsf_k3_tsne_ClusterMem.txt'.format(out_dir)
    # print(cluster_mem_file)
    predict_label = pd.read_csv(cluster_mem_file, delimiter='\t', index_col=0)
    # print(predict_label)
    merged = true_label.merge(predict_label, left_on='cell', right_on='ID')
    print(merged)
    ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
    print('{}'.format(ari))


if __name__ == "__main__":
    mica_data_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    # create_cmds()
    # summary()
    author = 'Buettner'
    calcu_ARI(author)

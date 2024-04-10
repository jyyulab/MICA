#!/usr/bin/env python3
import os
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score


out_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/outputs/SilverStd/PBMC_20k/activity/TALL_net'

summary_file = '{}/summary_ge.txt'.format(out_dir)
if os.path.isfile(summary_file):
    os.remove(summary_file)


true_label_file = '/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/PBMC_20k/' \
                  'PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
best_ari = 0.0


num_clusters = 10
for res in os.listdir(out_dir):
    if res.find('.txt') != -1:
        if res.find('clustering_UMAP_') != -1:
            print(res)
            with open(summary_file, 'w') as fout:
                predict_label = pd.read_csv('{}/{}'.format(out_dir, res), delimiter='\t', index_col=0)
                predict_num_clusters = len(set(predict_label['label']))
                if predict_num_clusters != num_clusters:
                    continue
                merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
                print('ARI: {}'.format(ari))
                # if ari > best_ari:
                #    best_ari = ari
                #   best_ari_str = '{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ari)
                # fout.write('{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ari))

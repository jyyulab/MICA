#!/usr/bin/env python3

from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd


true_label_file = './tests/answerkey/clustering/Buettner_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', dtype=int, header=0)
predict_label_file = './test_data/outputs/Buettner/clustering_UMAP_euclidean_20_3.txt'
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print('Buettner dataset ARI: {}'.format(ari))  # Expected ARI: 0.2662385902792365


true_label_file = './tests/answerkey/clustering/Yan_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', dtype=int, header=0)
predict_label_file = './test_data/outputs/Yan/Yan_k8_umap_ClusterMem.txt'
predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
merged = true_label.merge(predict_label, left_on='cell', right_index=True)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print('Yan dataset ARI: {}'.format(ari))    # Expected ARI: 0.9390370609566371

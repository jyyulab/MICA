#!/usr/bin/env python3

from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd

true_label_file = './tests/answerkey/clustering/kolod_true_label_1.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', dtype=int, header=0)
true_label = true_label.iloc[0]

predict_label_file = './test_data/outputs/clustering_umap_euclidean_12.txt'
predict_label = pd.read_csv(predict_label_file, delimiter='\t')
predict_label = predict_label['label']

ari = adjusted_rand_score(true_label, predict_label)
print('Kolod dataset ARI: {}'.format(ari))

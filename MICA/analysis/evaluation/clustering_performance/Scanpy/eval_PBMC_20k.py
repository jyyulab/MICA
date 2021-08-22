#!/usr/bin/env python3
from MICA.lib import preprocessing
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'PBMC_20k'
input_file = '{}/{}/{}/{}_MICA_input.h5ad'.format(root_dir, level, author, author)
num_clusters = 20

#%%
adata = preprocessing.read_preprocessed_mat(input_file)

#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#%%
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

#%%
sc.tl.pca(adata, svd_solver='arpack')

#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#%%
sc.tl.leiden(adata, resolution=0.08)
print(adata.obs['leiden'])

#%%
# sc.pl.umap(adata)

#%%
true_label_file = '{}/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
predict_label = adata.obs['leiden'].astype(int)
# predict_label.index = predict_label.index.astype(int)
predict_label.to_csv('{}/{}/{}/{}_predict_label.txt'.format(root_dir, level, author, author), sep='\t')

#%%
merged = true_label.merge(predict_label, left_on='cell', right_index=True)

#%%
ari = adjusted_rand_score(merged['label'], merged['leiden'])
print(ari)
#!/usr/bin/env python3
from MICA.lib import preprocessing
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'Chung'
input_file = '{}/{}/{}/GSE75688_{}_preprocessed_317.h5ad'.format(root_dir, level, author, author)
num_clusters = 4

#%%
adata = preprocessing.read_preprocessed_mat(input_file)

#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
true_label = pd.read_csv(true_label_file, delimiter=',', header=0)

#%%
adata = adata[list(true_label['cell_name']), :]

#%%
results_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE75688_Chung/' \
               'GSE75688_Chung_preprocessed_317.h5ad'
adata.write(results_file)


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
sc.tl.leiden(adata, resolution=0.3)
print(adata.obs['leiden'])


#%%
sc.tl.umap(adata)

#%%
sc.pl.umap(adata)

#%% UMAP scatter plot
df_umap = pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'])

#%%
df_umap['label'] = list(adata.obs['leiden'].astype(int))

#%%
df_umap.index = adata.obs.index

#%%
df_umap.to_csv('/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Chung/Scanpy/Chung_Scanpy_UMAP.txt',
               sep='\t')

#%%
# true_label_file = '{}/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
true_label = pd.read_csv(true_label_file, delimiter=',', header=0)

#%%
predict_label = adata.obs['leiden'].astype(int)
# predict_label.index = predict_label.index.astype(int)

#%%
merged = true_label.merge(predict_label, left_on='cell_name', right_index=True)

#%%
ari = adjusted_rand_score(merged['cell_type'], merged['leiden'])
ari

#%%
ami = adjusted_mutual_info_score(merged['cell_type'], merged['leiden'])
ami

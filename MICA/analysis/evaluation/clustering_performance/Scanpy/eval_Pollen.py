#!/usr/bin/env python3
from MICA.lib import preprocessing
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'GoldenStd'
author = 'Pollen'
input_file = '{}/{}/{}/{}_MICA_input_246.h5ad'.format(root_dir, level, author, author)
num_clusters = 11

#%%
adata = preprocessing.read_preprocessed_mat(input_file)

#%%
sc.pl.highest_expr_genes(adata, n_top=20, )

#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

#%%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

#%%
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

#%%
sc.pp.log1p(adata)

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
sc.tl.leiden(adata, resolution=2.7)
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
df_umap.to_csv('/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Pollen/Scanpy/Pollen_Scanpy_UMAP.txt',
               sep='\t')

#%%
true_label_file = '{}/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
predict_label = adata.obs['leiden'].astype(int)
predict_label.index = predict_label.index

#%%
merged = true_label.merge(predict_label, left_on='cell', right_index=True)

#%%
ari = adjusted_rand_score(merged['label'], merged['leiden'])
ari

#%%
ami = adjusted_mutual_info_score(merged['label'], merged['leiden'])
ami

#%%
df = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
df.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_MICA_input_246.txt', sep='\t')

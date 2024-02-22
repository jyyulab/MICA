#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics.cluster import adjusted_rand_score


#%%
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%%
pbmc20k_h5ad = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
               'PBMC_sorted_20K_preprocessed.h5ad'
adata = ad.read_h5ad(filename=pbmc20k_h5ad)

#%%
sc.pl.highest_expr_genes(adata, n_top=20)

#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

#%%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

#%%
sc.pl.highly_variable_genes(adata)

#%%
adata.raw = adata

#%%
adata = adata[:, adata.var.highly_variable]

#%%
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%%
sc.pp.scale(adata, max_value=10)

#%%
sc.tl.pca(adata, svd_solver='arpack')

#%%
sc.pl.pca(adata, color='CST3')

#%%
sc.pl.pca_variance_ratio(adata, log=True)

#%%
pbmc20k_h5ad_high_variable = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                             'PBMC_sorted_20K_preprocessed_high_variable.h5ad'

#%%
adata.write(pbmc20k_h5ad_high_variable)

#%%
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)

#%%
sc.tl.umap(adata)

#%%
sc.tl.leiden(adata, resolution=0.14)

#%%
sc.pl.umap(adata, color=['leiden','CST3'])

#%%
adata.write(pbmc20k_h5ad_high_variable)

#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                  'PBMC_sorted_20K_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
cell_type_label_dict = dict()
for i, v in enumerate(set(true_label['type'])):
    cell_type_label_dict[v] = i
labels = [cell_type_label_dict[ct] for ct in true_label['type']]
true_label['label'] = labels

#%% Scanpy
predict_label = adata.obs.leiden
adjusted_rand_score(true_label['label'], predict_label)
# ARI: 0.08013603368345233


#%% Seurat
predict_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                     'seurat_predicted_label.csv'
predict_label = pd.read_csv(predict_label_file, delimiter=',', header=0)
adjusted_rand_score(true_label['label'], predict_label['label'])
# ARI: 0.5777771977018957

#%%
adata = ad.read_h5ad(filename=pbmc20k_h5ad_high_variable)

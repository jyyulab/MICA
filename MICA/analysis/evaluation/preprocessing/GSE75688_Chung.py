#!/usr/bin/env python3
"""
Preprocessing the raw data to create an AnnData object.
10X data source: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats as stats
import anndata
import logging
from collections import Counter
from sklearn.metrics.cluster import adjusted_rand_score


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
logging.basicConfig(level=logging.INFO)


#%%
raw_TPM_path = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE75688_Chung/raw/' \
               'GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_no_pool.txt'
raw_TPM_df = pd.read_csv(raw_TPM_path, delimiter='\t', index_col=0)


#%%
raw_TPM_df_T = raw_TPM_df.T


#%%
raw_TPM_df_T.to_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE75688_Chung/raw/' \
                    'GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_no_pool_T.txt')


#%%
adata = sc.read_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE75688_Chung/raw/' \
                    'GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix_no_pool_T.txt', first_column_names=True)


#%% Optional, see Scanpy clustering tutorial for details
adata.var_names_make_unique()


#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
true_label = pd.read_csv(true_label_file, delimiter=',', header=0)

#%%
adata = adata[list(true_label['cell_name']), :]


#%%
sc.pl.highest_expr_genes(adata, n_top=20, )

#%%
# sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%%
sc.pl.highest_expr_genes(adata, n_top=20, )

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
sc.pp.log1p(adata)

#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

#%%
sc.pl.highly_variable_genes(adata)

#%%
adata.raw = adata

#%%
adata = adata[:, adata.var.highly_variable]

#%%
sc.pp.scale(adata, max_value=10)


#%%
sc.tl.pca(adata, svd_solver='arpack')


#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#%%
sc.tl.leiden(adata, resolution=0.1)
print(adata.obs['leiden'])

#%%
sc.tl.umap(adata)

#%%
sc.pl.umap(adata)




#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/scGNN/Chung_cell_label.csv'
true_label = pd.read_csv(true_label_file, delimiter=',', header=0)


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'Chung'
num_clusters = 4

#%%
predict_label = adata.obs['leiden'].astype(int)
# predict_label.index = predict_label.index.astype(int)
predict_label.to_csv('{}/{}/{}/{}_predict_label.txt'.format(root_dir, level, author, author), sep='\t')


#%%
merged = true_label.merge(predict_label, left_on='cell_name', right_index=True)

#%%
ari = adjusted_rand_score(merged['cell_type'], merged['leiden'])
print(ari)

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
sc.pl.highest_expr_genes(adata, n_top=20, )

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
true_label_df = pd.read_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/'
                            'GSE75688_Chung/GSE75688_true_label.txt', sep='\t')

#%%
true_label_df['sample']

#%%
adata = adata[true_label_df['sample'], adata.var.highly_variable]

#%%
results_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE75688_Chung/' \
               'GSE75688_Chung_preprocessed.h5ad'
adata.write(results_file)

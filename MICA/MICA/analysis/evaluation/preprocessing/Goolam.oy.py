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
raw_count_path = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Goolam/' \
                 'Goolam_et_al_2015_count_table.tsv'
raw_count_df = pd.read_csv(raw_count_path, delimiter='\t', index_col=0)

#%%
raw_count_df_T = raw_count_df.T

#%%
raw_count_df_T.to_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Goolam/' \
                      'Goolam_et_al_2015_count_table_T.tsv')

#%%
adata = sc.read_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Goolam/' \
                    'Goolam_et_al_2015_count_table_T.tsv', first_column_names=True)

#%% Optional, see Scanpy clustering tutorial for details
adata.var_names_make_unique()

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

#%% Only 124 cells, no need to filter cells
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

#%%
sc.pp.normalize_total(adata, target_sum=1e4)

#%%
sc.pp.log1p(adata)

#%%
results_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Goolam/Goolam_MICA_input.h5ad'
adata.write(results_file)

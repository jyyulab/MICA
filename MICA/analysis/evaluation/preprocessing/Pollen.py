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
from MICA.lib import preprocessing as pp

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
logging.basicConfig(level=logging.INFO)

#%%
adata = pp.read_preprocessed_mat('/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_MICA_input.txt')

#%%
adata.var_names_make_unique()

#%%
sc.pp.filter_genes(adata, min_cells=10)

#%%
sc.pp.filter_cells(adata, min_genes=6000.000)
#%%
# sc.pp.filter_cells(adata, max_genes=10000.000)

#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

#%%
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

#%%
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]

#%%
sc.pp.log1p(adata)

#%%
adata.obs.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_obs.txt', sep='\t')

#%%
adata.write_h5ad(filename='/Users/lding/Documents/MICA/Datasets/HPC/GoldenStd/Pollen/Pollen_MICA_input_246.h5ad')

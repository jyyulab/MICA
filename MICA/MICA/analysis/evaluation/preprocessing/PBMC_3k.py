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
adata = sc.read_10x_mtx(
    'test_data/inputs/10x/PBMC/3k/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

#%% Optional, see Scanpy clustering tutorial for details
adata.var_names_make_unique()


#%%
sc.pl.highest_expr_genes(adata, n_top=20, )


#%% Compute QC metrics
# annotate the group of mitochondrial genes as 'mt' or 'MT'
adata.var['mt'] = np.logical_or(adata.var_names.str.startswith('mt-'), adata.var_names.str.startswith('MT-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


#%% Set total UMI cutoff to be max of median +/- 3 * MAD
umi_cutoff_low = max([np.median(adata.obs['total_counts']) -
                      3 * stats.median_absolute_deviation(adata.obs['total_counts']), 100])
umi_cutoff_high = np.median(adata.obs['total_counts']) + 3 * stats.median_absolute_deviation(adata.obs['total_counts'])
logging.info('Total UMI per cell low/high cutoffs: {}/{}'.format(umi_cutoff_low, umi_cutoff_high))


#%% Filter cells based on total UMIs per cell
sc.pp.filter_cells(adata, min_counts=umi_cutoff_low)
sc.pp.filter_cells(adata, max_counts=umi_cutoff_high)


#%% Set number of genes per cell to be max of median +/- 3 * MAD
min_gene_cutoff = max([np.floor(np.median(adata.obs['n_genes_by_counts']) -
                      3 * stats.median_absolute_deviation(adata.obs['n_genes_by_counts'])), 50])
max_gene_cutoff = min([np.floor(np.median(adata.obs['n_genes_by_counts']) +
                      3 * stats.median_absolute_deviation(adata.obs['n_genes_by_counts'])), 2500])
logging.info('Minimum/Maximum genes per cell cutoffs: {}/{}'.format(min_gene_cutoff, max_gene_cutoff))

#%% Filter cells based on total genes expressed in a cell
sc.pp.filter_cells(adata, min_genes=min_gene_cutoff)
sc.pp.filter_cells(adata, max_genes=max_gene_cutoff)


#%% Filter genes expressed less than certain percentage of cells
percent_cutoff = 0.005        # default cutoff is 0.005, i.e., expressing less than 0.5% of cells will be filtered
num_cells = int(adata.shape[0] * percent_cutoff)
sc.pp.filter_genes(adata, min_cells=num_cells)

#%% Filter cells with >5% mitochondria counts
adata = adata[adata.obs.pct_counts_mt < 5, :]

#%% Normalize the count matrix to 1,000,000 per cell and Logarithmize the data
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata, base=2)

#%% Save a copy of pre-processed data
adata.raw = adata


#%% The following 3 steps are suggested by Seurat and Scanpy. Their performance on MICA is under evaluation.
# (Optional) Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%% (Optional) Scale the data to unit variance.
# sc.pp.scale(adata, max_value=10)

#%% (Optional) Identify and plot highly variable genes
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc.pl.highly_variable_genes(adata)
# Actually do the filtering on the matrix
# adata = adata[:, adata.var.highly_variable]


#%% Save the pre-processed sparse matrix
preprocessed_results = './test_data/inputs/10x/PBMC/3k/pre-processed/pbmc3k_preprocessed.h5ad'
adata.write(preprocessed_results)

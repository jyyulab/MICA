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
h5_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/LabelFree/' \
          '10X_1.3_Million_Brain_Cells_from_E18_Mice/1M_neurons_filtered_gene_bc_matrices_h5.h5'
adata = sc.read_10x_h5(h5_path, genome='mm10')

#%%
adata.var_names_make_unique()

#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%%
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4)
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, save='n_genes_by_counts.pdf')
sc.pl.violin(adata, ['total_counts'], jitter=0.4, save='total_counts.pdf')
sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, save='pct_counts_mt.pdf')

#%%
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

#%%
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

#%%
sc.pl.highly_variable_genes(adata, save='highly_variable_genes.pdf')

#%%
results_file = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/LabelFree/' \
               '10X_1.3_Million_Brain_Cells_from_E18_Mice/1M_neurons_filtered_gene_bc_matrices.h5ad'
adata.write(results_file)

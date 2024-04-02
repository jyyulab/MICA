#!/usr/bin/env python3
"""
Preprocessing the raw data to create an AnnData object.
10X data source: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc33k
"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.stats as stats
import anndata
import logging
from collections import Counter
import memory_profiler

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
logging.basicConfig(level=logging.INFO)


#%%
adata = sc.read_h5ad('/Users/lding/Documents/MICA/Datasets/filtered_gene_bc_matrices/hg19/pbmc33k_preprocessed.h5ad')

#%%
import time
@profile
def mem_profile():
    start = time.time()
    sc.tl.pca(adata, svd_solver='arpack')
    end = time.time()
    print(end - start)

    #%%
    sc.pl.pca(adata, color='CST3')

    #%%
    start = time.time()
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    end = time.time()
    print(end - start)

    #%%
    start = time.time()
    sc.tl.leiden(adata)
    end = time.time()
    print(end - start)


if __name__ == '__main__':
    mem_profile()

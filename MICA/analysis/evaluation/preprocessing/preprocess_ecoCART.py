#!/usr/bin/env python3
"""
Preprocessing the raw data to create an AnnData object.
10X data source: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc33k
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
raw_count_path = '/Users/lding/Desktop/ecoCART_double_construct_bothdonor_MICAfilt_except_week19.txt'
adata = sc.read_csv(raw_count_path, first_column_names=True, delimiter='\t')

#%%
results_file = '/Users/lding/Desktop/ecoCART_double_construct_bothdonor_MICAfilt_except_week19.h5ad'
adata.write(results_file)

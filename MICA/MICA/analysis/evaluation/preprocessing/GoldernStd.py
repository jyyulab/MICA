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
df = pd.read_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/buettner/Buttner_MICA_input.txt',
                 index_col=0, delimiter='\t', header=0)

#%%
adata = anndata.AnnData(X=df)


#%%
results_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/buettner/Buttner_MICA_input.h5ad'
adata.write(results_file)

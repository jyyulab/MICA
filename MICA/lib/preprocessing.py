#!/usr/bin/env python3

import pandas as pd
import anndata


def read_preprocessed_mat(in_file):
    """ Read in preprocessed matrix file (h5ad file or tab-delimited text file) into a dataframe."""
    if in_file.endswith('.txt'):
        frame = pd.read_csv(in_file, sep="\t", index_col=0).iloc[:, 0:]
    if in_file.endswith('.h5ad') or in_file.endswith('.h5'):
        adata = anndata.read_h5ad(in_file)
        frame = adata.to_df()
    return frame

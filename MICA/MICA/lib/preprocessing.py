#!/usr/bin/env python3

import sys
import pandas as pd
import anndata
import numpy as np
from natsort import natsorted


def read_preprocessed_mat(in_file, start_row=1):
    """ Read in preprocessed matrix file (h5ad file or tab-delimited text file) into a dataframe."""
    if in_file.endswith('.txt'):
        adata = anndata.read_text(in_file, delimiter='\t', first_column_names=True)
        if adata.obs.index[0] == 'ID' or adata.var.index[0] == 0:
            adata = adata[start_row:, :]    # 1st row used as header by default
    elif in_file.endswith('.csv'):
        adata = anndata.read_csv(in_file, first_column_names=True)
        if adata.obs.index[0] == 'ID' or adata.var.index[0] == 0:
            adata = adata[start_row:, :]    # 1st row used as header by default
    elif in_file.endswith('.h5ad') or in_file.endswith('.h5'):
        adata = anndata.read_h5ad(in_file)
    else:
        sys.exit('Error - invalid input file format: {}'.format(in_file))
    return adata


def read_write_mat(in_file, out_file_name, start_row=1):
    """ Read in preprocessed matrix file (h5ad file or tab-delimited text file) into a dataframe
    and write the entire matrix as slice_0 for further slicing.
    """
    adata = read_preprocessed_mat(in_file, start_row)
    frame = adata.to_df()
    frame.to_hdf(out_file_name + ".h5.tmp", "slice_0")


def write_h5(adata, out_h5_file):
    """ Output annData as a h5ad file in the given output dir. """
    adata.write(out_h5_file)


def write_h5_clustering(adata, out_h5_file, partition, resolution):
    """ Output annData with clustering information as a h5ad file in the given output dir. """
    labels = [x + 1 for x in partition.values()]
    adata.obs['louvain'] = pd.Categorical(
        values=labels,
        categories=natsorted(map(str, np.unique(labels))),
    )
    # store information on the clustering parameters
    adata.uns['louvain'] = {}
    adata.uns['louvain']['params'] = dict(
        resolution=resolution,
    )
    adata.write(out_h5_file)
    return adata

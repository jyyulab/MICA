#!/usr/bin/env python3

import sys
import pandas as pd
import anndata
import numpy as np
from natsort import natsorted
from anndata import AnnData
from typing import Union, Optional, Mapping
from typing import Iterable, Iterator, Generator
from os import PathLike, fspath


def read_preprocessed_mat(in_file):
    """ Read in preprocessed matrix file (h5ad file or tab-delimited text file) into a dataframe."""
    if in_file.endswith('.txt'):
        adata = anndata.read_csv(in_file, delimiter='\t', first_column_names=True)
        adata = adata[1:, :]    # 1st row used as header
    elif in_file.endswith('.csv'):
        adata = anndata.read_csv(in_file)
    elif in_file.endswith('.h5ad') or in_file.endswith('.h5'):
        adata = anndata.read_h5ad(in_file)
    else:
        sys.exit('Error - invalid input file format: {}'.format(in_file))
    return adata


def write_h5(adata, out_h5_file, partition, resolution):
    """ Output annData as a h5ad file in the given output dir. """
    labels = [x + 1 for x in partition.values()]
    adata.obs['louvain'] = pd.Categorical(
        values=labels,
        categories=natsorted(map(str, np.unique(labels))),
    )
    # store information on the clustering parameters
    adata.uns['leiden'] = {}
    adata.uns['leiden']['params'] = dict(
        resolution=resolution,
    )
    adata.write(out_h5_file)
    return adata

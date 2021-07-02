#!/usr/bin/env python3

from MICA.lib import preprocessing
from collections import Counter
import pandas as pd
import h5py
import anndata


#%%
infile = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/PBMC_sorted_20K_preprocessed.h5ad'
adata3 = preprocessing.read_preprocessed_mat(infile)

#%%
print(adata3.X.shape)
print(adata3.obs.shape)

#%%
infile = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE71585_Tasic/GSE71585_MICA_input.txt'
adata1 = anndata.read_csv(infile, delimiter='\t', first_column_names=True)

#%%
print(adata1.X.shape)
print(adata1.obs.shape)

#%%
start_row = 1
adata = adata1[start_row:, :]

#%%
infile = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Usoskin/Usoskin_MICA_input.txt'
adata2 = anndata.read_csv(infile, delimiter='\t', first_column_names=True)

#%%
print(adata2.X.shape)
print(adata2.obs.shape)

#%%
start_row = 1
adata = adata[start_row:, :]

#%%
frame = adata.to_df()

#%%
f = h5py.File(infile)

#%%
print(f.keys())

#%%
mat = f['input']

#%%
adata = preprocessing.read_preprocessed_mat(infile)

#%%
frame.to_hdf("/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Klein/Klein.h5.tmp", "slice_0")

#%%
freq = Counter(adata.var.index)

#%%
for key in freq:
    if freq[key] > 1:
        print(key)

#%%
preprocessing.read_write_mat(infile, 'Klein')

#%%
infile = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Klein/Klein.txt'
frame = pd.read_csv(infile, delimiter='\t')

#%%
frame_transpose = frame.T

#%%
frame_transpose.to_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/Klein/Klein_MICA_input2.txt',
                       sep='\t')

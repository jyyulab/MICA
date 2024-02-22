#!/usr/bin/env python3

import os
import pandas as pd

#%%
ziesel_true_label_path = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE60361_Ziesel/' \
                         'Ziesel_true_label.txt'
ziesel_tl_df = pd.read_csv(ziesel_true_label_path, delimiter='\t')


#%%
ziesel_tl_df_T = ziesel_tl_df.T

#%%
ziesel_tl_df_T.to_csv('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/GSE60361_Ziesel/'
                      'Ziesel_true_label_T.txt', sep='\t')

#%%
import scanpy as sc

#%%
adata = sc.read_h5ad('/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/'
                     'PBMC_sorted_20K_preprocessed.h5ad')

#%%
index_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/index.txt'
with open(index_file, 'w') as fout:
    for barcode in list(adata.obs.index):
        fout.write('{}\n'.format(barcode))

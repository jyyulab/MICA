#!/usr/bin/env python3
import scanpy as sc
import numpy as np
import anndata as ad
from collections import Counter


#%%
adata = sc.read_h5ad("/Users/lding/Documents/MICA/Datasets/revision/covid/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection.h5ad")

#%%
print(adata)

#%%
adata.obs

#%%
Counter(adata.X[0])

#%%
adata.X = adata.layers['matrix']

#%%
Counter(adata.X[0])

#%%
adata.X = adata.X.astype(int)

#%%
Counter(adata.X[0])

#%%
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

#%%
adata.write_h5ad("/Users/lding/Documents/MICA/Datasets/revision/covid/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection_CPM_log1p.h5ad",
                 compression="gzip")

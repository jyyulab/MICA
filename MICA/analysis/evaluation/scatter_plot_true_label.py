#!/usr/bin/env python3
import pandas as pd
import copy
from MICA.lib import visualize as vs

#%%
mica_dir = '/Users/lding/Documents/MICA/outputs'
ge_umap_out_file = '{}/PBMC_20k_default_parameters/40_1.4/clustering_umap_euclidean_40_1.4.txt'.format(mica_dir)
mds_umap_out_file = '{}/PBMC_20k_MDS/cwl_lsf_k12_tsne_ClusterMem.txt'.format(mica_dir)
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                  'PBMC_sorted_20K_true_label.txt'

#%%
ge_umap = pd.read_csv(ge_umap_out_file, delimiter='\t', index_col=0)
ge_umap_true_label = copy.deepcopy(ge_umap)

#%%
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_umap_true_label['label'] = true_label

#%%
out_png_file = '{}/PBMC_20k_default_parameters/40_1.4/clustering_umap_euclidean_40_1.4_true_label.png'.format(mica_dir)
vs.scatter_plot(ge_umap_true_label, out_png_file)


#%%
mds_umap = pd.read_csv(mds_umap_out_file, delimiter='\t', index_col=0)
mds_umap_true_label = copy.deepcopy(mds_umap)
mds_umap_true_label['label'] = true_label

#%%
out_png_file = '{}/PBMC_20k_MDS/cwl_lsf_k12_tsne_ClusterMem_true_label.png'.format(mica_dir)
vs.scatter_plot(mds_umap_true_label, out_png_file)

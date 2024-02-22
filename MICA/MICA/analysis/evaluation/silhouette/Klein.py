#!/usr/bin/env python3

from MICA.lib import visualize as vi
import pandas as pd


#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Klein/Klein_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
# true_label['label'] = true_label['label']+1


#%%
umap_embed_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/mds_1_3_silhouette/clustering_umap_euclidean_19.txt'
umap_embed = pd.read_csv(umap_embed_file, sep='\t', index_col=0)
# umap_embed['label'] = umap_embed['label'] + 1


#%%
merged = true_label.merge(umap_embed, left_on='cell', right_index=True)


#%%
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/mds_1_3_silhouette'
labels = merged['label_x']


#%%
label_dict = dict()
for i, l in enumerate(set(labels)):
    label_dict[l] = i + 1

#%%
num_labels = []
for l in labels:
    num_labels.append(label_dict[l])

#%%
merged['num_labels'] = num_labels

#%%
labels = merged['num_labels']


#%%
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(merged.loc[:, ['X', 'Y']], labels)
print(silhouette_avg)
# sample_silhouette_values = silhouette_samples(merged.loc[:, ['X', 'Y']], labels)

#%%
vi.silhouette_plot(labels, merged.loc[:, ['X', 'Y']], 4, silhouette_avg, out_dir, ss_lower_bound=-0.6,
                   ss_upper_bound=1.0)





#%% Seurat
umap_embed_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Seurat/Klein_Seurat_UMAP.txt'
umap_embed = pd.read_csv(umap_embed_file, sep='\t', index_col=0)

#%%
umap_embed.index = list(range(0, 2717))

#%%
merged = true_label.merge(umap_embed, left_index=True, right_index=True)

#%%
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Seurat'

#%%
labels = merged['label']

#%%
label_dict = dict()
for i, l in enumerate(set(labels)):
    label_dict[l] = i + 1

#%%
num_labels = []
for l in labels:
    num_labels.append(label_dict[l])

#%%
merged['num_labels'] = num_labels

#%%
labels = merged['num_labels']

#%%
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(merged.loc[:, ['UMAP_1', 'UMAP_2']], labels)
print(silhouette_avg)

#%%
vi.silhouette_plot(labels, merged.loc[:, ['UMAP_1', 'UMAP_2']], 4, silhouette_avg, out_dir, ss_lower_bound=-0.6,
                   ss_upper_bound=1.0)



#%% Scanpy
umap_embed_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy/Klein_Scanpy_UMAP.txt'
umap_embed = pd.read_csv(umap_embed_file, sep='\t', index_col=0)

#%%
merged = true_label.merge(umap_embed, left_index=True, right_index=True)

#%%
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy'

#%%
labels = merged['label_x']

#%%
label_dict = dict()
for i, l in enumerate(set(labels)):
    label_dict[l] = i + 1

#%%
num_labels = []
for l in labels:
    num_labels.append(label_dict[l])

#%%
merged['num_labels'] = num_labels

#%%
labels = merged['num_labels']

#%%
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(merged.loc[:, ['X', 'Y']], labels)
print(silhouette_avg)


#%%
vi.silhouette_plot(labels, merged.loc[:, ['X', 'Y']], 4, silhouette_avg, out_dir, ss_lower_bound=-0.6,
                   ss_upper_bound=1.0)





#%% SC3
umap_embed_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/SC3/Klein_SC3_euclidean_laplacian_UMAP.txt'
umap_embed = pd.read_csv(umap_embed_file, sep='\t', index_col=0)

#%%
umap_embed.index -= 1

#%%
merged = true_label.merge(umap_embed, left_index=True, right_index=True)

#%%
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/SC3'

#%%
labels = merged['label']

#%%
label_dict = dict()
for i, l in enumerate(set(labels)):
    label_dict[l] = i + 1

#%%
num_labels = []
for l in labels:
    num_labels.append(label_dict[l])

#%%
merged['num_labels'] = num_labels

#%%
labels = merged['num_labels']

#%%
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(merged.loc[:, ['V1', 'V2']], labels)
print(silhouette_avg)


#%%
vi.silhouette_plot(labels, merged.loc[:, ['V1', 'V2']], 4, silhouette_avg, out_dir, ss_lower_bound=-0.6,
                   ss_upper_bound=1.0)



#%%
out_pdf_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/SC3/' \
               'Klein_SC3_euclidean_laplacian_UMAP.pdf'
umap_embed['label'] = labels

#%%
umap_embed = umap_embed.rename(columns={'V1': 'X', 'V2': 'Y'})

#%%
vi.scatter_plot(umap_embed, out_pdf_file, method='umap')

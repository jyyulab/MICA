#!/usr/bin/env python3
from MICA.lib import preprocessing
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'Klein'
input_file = '{}/{}/{}/{}_MICA_input.txt'.format(root_dir, level, author, author)
num_clusters = 4

#%%
adata = preprocessing.read_preprocessed_mat(input_file)

#%%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#%%
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

#%%
sc.tl.pca(adata, svd_solver='arpack')

#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#%%
sc.tl.leiden(adata, resolution=0.1)
print(adata.obs['leiden'])

#%%
sc.tl.umap(adata)

#%%
sc.pl.umap(adata)

#%% UMAP scatter plot
df_umap = pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'])

#%%
df_umap['label'] = list(adata.obs['leiden'].astype(int))

#%%
df_umap.to_csv('/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy/Klein_Scanpy_UMAP.txt',
               sep='\t')



#%%
true_label_file = '{}/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
predict_label = adata.obs['leiden'].astype(int)
# predict_label.index = predict_label.index.astype(int)

#%%
merged = true_label.merge(predict_label, left_on='cell', right_index=True)

#%%
ari = adjusted_rand_score(merged['label'], merged['leiden'])
ari

#%%
ami = adjusted_mutual_info_score(merged['label'], merged['leiden'])
ami

#%%
sc.tl.umap(adata)

#%%
sc.pl.umap(adata)



#%% UMAP scatter plot
df_umap = pd.DataFrame(adata.obsm['X_umap'], columns=['X', 'Y'])

#%%
df_umap['label'] = list(adata.obs['leiden'].astype(int))
df_umap['label'] = df_umap['label'] + 1

#%%
df_umap.to_csv('/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy/Klein_Scanpy_UMAP.txt',
               sep='\t')

#%% UMAP scatter plot
from MICA.lib import visualize as vi
out_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy/Klein_Scanpy_UMAP.pdf'
vi.scatter_plot(df_umap, out_file, marker_size=1.0, marker="o", method='UMAP', marker_scale=10.0)


#%% silhouette score
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
print(silhouette_avg)
# 0.25015703


#%%
df_pca = pd.DataFrame(adata.obsm['X_pca'])
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Klein/Scanpy'
labels = adata.obs['leiden']
labels = labels.astype(int)

#%%
df_pca = df_pca.set_index(labels.index)

#%%
from sklearn.metrics import silhouette_samples
sample_silhouette_values = silhouette_samples(df_pca, labels)
from MICA.lib import visualize as vi
vi.silhouette_plot(labels, df_pca.to_numpy(), 4, 0.25, out_dir, ss_lower_bound=-0.6, ss_upper_bound=0.8)

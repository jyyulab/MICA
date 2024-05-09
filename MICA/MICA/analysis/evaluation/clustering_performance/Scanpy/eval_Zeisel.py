#!/usr/bin/env python3
from MICA.lib import preprocessing
import scanpy as sc
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'Zeisel'
num_clusters = 7


#%% Read raw data
input_file = '{}/{}/{}/GSE60361_C1-3005-Expression.txt'.format(root_dir, level, author)
adata = sc.read_csv(input_file, first_column_names=True, delimiter='\t')

#%%
# input_file = '{}/{}/{}/Zeisel_preprocessed.h5ad'.format(root_dir, level, author)
# adata = preprocessing.read_preprocessed_mat(input_file)


#%%
adata = adata.transpose()


#%%
adata.var_names_make_unique()


#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


#%%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.3, multi_panel=True)


#%%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


#%%
# adata2 = adata[adata.obs.n_genes_by_counts < 8000, :]

#%%
adata = adata[adata.obs.pct_counts_mt < 5, :]


#%%
sc.pp.normalize_total(adata, target_sum=1e4)

#%%
sc.pp.log1p(adata)

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
sc.tl.leiden(adata, resolution=0.125)
print(adata.obs['leiden'])

#%%
sc.tl.umap(adata)

#%%
sc.pl.umap(adata)



#%%
input_file = '{}/{}/{}/{}_MICA_input.txt'.format(root_dir, level, author, author)
num_clusters = 7

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

#%%
true_label_file = '{}/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
predict_label = adata.obs['leiden'].astype(int)
# predict_label.index = predict_label.index.astype(int)

#%%
true_label['cell'] = true_label['cell'].str.lstrip('X')

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

#%%
df_umap.to_csv('/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Zeisel/Scanpy/Zeisel_Scanpy_UMAP.txt',
               sep='\t')

#%% UMAP scatter plot
from MICA.lib import visualize as vi
df_umap['label'] += 1
out_file = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Zeisel/Scanpy/Zeisel_Scanpy_UMAP.pdf'
vi.scatter_plot(df_umap, out_file, marker_size=1.0, marker="o", method='UMAP', marker_scale=10.0)


#%% silhouette score
from sklearn.metrics import silhouette_score
silhouette_avg = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
print(silhouette_avg)




#%%
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def silhouette_plot(labels, frame_dr, num_clusters, silhouette_avg, out_dir, ss_lower_bound=-0.1, ss_upper_bound=1.0):
    """ Draw a silhouette plot.
    Args:
        labels (array): array-like clustering results, {sample index: cluster label}
        frame_dr (ndarray): matrix after dimension reduction
        num_clusters (int): number of clusters in labels
        silhouette_avg (float): silhouette score
        out_dir (dir): path to output folder
        ss_lower_bound (float): The silhouette coefficient can range from -1, 1. This parameter sets the lower
        limit for plotting.
        ss_upper_bound (float): The silhouette coefficient can range from -1, 1. This parameter sets the upper
        limit for plotting.
    Returns:
        PDF image of silhouette plot
    """
    # Compute the silhouette scores for each cell
    sample_silhouette_values = silhouette_samples(frame_dr, labels)
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 7)
    y_lower = 10
    if min(labels) == 0:
        labels = labels + 1
    for i in range(1, num_clusters + 1):
        # Aggregate the silhouette scores for cells belonging to cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[labels == i]
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        # The silhouette coefficient can range from -1, 1
        ax.set_xlim([ss_lower_bound, ss_upper_bound])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax.set_ylim([0, len(frame_dr) + (num_clusters + 1) * 10])

        color = cm.nipy_spectral(float(i) / num_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                         0, ith_cluster_silhouette_values,
                         facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks(np.arange(ss_lower_bound, ss_upper_bound+0.1, 0.2))
    out_png_file = '{}/silhouette_{}_{}.pdf'.format(out_dir, frame_dr.shape[1], num_clusters)
    plt.savefig(out_png_file, bbox_inches="tight")
    return sample_silhouette_values



#%%
df_pca = pd.DataFrame(adata.obsm['X_pca'])
out_dir = '/Users/lding/Documents/MICA/Manuscript/Figures/Silhouette/Zeisel/Scanpy'
labels = adata.obs['leiden']
labels = labels.astype(int)

#%%
df_pca = df_pca.set_index(labels.index)

#%%
from sklearn.metrics import silhouette_samples
sample_silhouette_values = silhouette_samples(df_pca, labels)
silhouette_plot(labels, df_pca.to_numpy(), 7, 0.0106, out_dir, ss_lower_bound=-0.6, ss_upper_bound=0.8)

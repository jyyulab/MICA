#!/usr/bin/env python3
import pandas as pd
import copy
from MICA.lib import visualize as vs



#%% Chung
mica_dir = '/Users/lding/Documents/MICA/outputs'
ge_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/clustering_UMAP_euclidean_60_2.2.txt'
mds_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/Chung_k5_umap_ClusterMem.txt'
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Chung/Chung_true_label.txt'
seurat_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/Chung_k5_seurat_umap.csv'

#%%
ge_umap = pd.read_csv(ge_umap_out_file, delimiter='\t', index_col=0)
ge_umap_true_label = copy.deepcopy(ge_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_umap_true_label['label'] = true_label

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/clustering_UMAP_euclidean_60_2.2_true_label.png'
vs.scatter_plot(ge_umap_true_label, out_png_file, method='umap', marker_size=2.0)

#%%
mds_umap = pd.read_csv(mds_umap_out_file, delimiter='\t', index_col=0)
mds_umap_true_label = copy.deepcopy(mds_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
mds_umap_true_label['label'] = true_label

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/Chung_k5_umap_ClusterMem_true_label.png'
vs.scatter_plot(mds_umap_true_label, out_png_file, method='umap', marker_size=2.0)

#%%
seurat_umap = pd.read_csv(seurat_umap_out_file, delimiter=',', index_col=0)
seurat_umap_true_label = copy.deepcopy(seurat_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
seurat_umap_true_label['label'] = true_label

#%%
seurat_umap_true_label = seurat_umap_true_label.rename(columns={'UMAP_1': 'X', 'UMAP_2': 'Y'})
seurat_umap_true_label = seurat_umap_true_label.loc[~seurat_umap_true_label['label'].isna()]

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE75688_Chung/Chung_k5_seurat_umap_true_label.png'
vs.scatter_plot(seurat_umap_true_label, out_png_file, method='umap', marker_size=2.0)





#%% Tasic
mica_dir = '/Users/lding/Documents/MICA/outputs'
ge_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/clustering_UMAP_euclidean_20_2.2.txt'
mds_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/Tasic_k8_umap_ClusterMem.txt'
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/Tasic/Tasic_true_label.txt'
seurat_umap_out_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/Tasic_k8_seurat_umap.csv'

#%%
ge_umap = pd.read_csv(ge_umap_out_file, delimiter='\t', index_col=0)
ge_umap_true_label = copy.deepcopy(ge_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_umap_true_label['label'] = true_label

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/clustering_UMAP_euclidean_20_2.2_true_label.png'
vs.scatter_plot(ge_umap_true_label, out_png_file, method='umap', marker_size=2.0)

#%%
mds_umap = pd.read_csv(mds_umap_out_file, delimiter='\t', index_col=0)
mds_umap_true_label = copy.deepcopy(mds_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
mds_umap_true_label['label'] = true_label

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/Tasic_k8_umap_ClusterMem_true_label.png'
vs.scatter_plot(mds_umap_true_label, out_png_file, method='umap', marker_size=2.0)

#%%
seurat_umap = pd.read_csv(seurat_umap_out_file, delimiter=',', index_col=0)
seurat_umap_true_label = copy.deepcopy(seurat_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
seurat_umap_true_label['label'] = true_label

#%%
seurat_umap_true_label = seurat_umap_true_label.rename(columns={'UMAP_1': 'X', 'UMAP_2': 'Y'})

#%%
out_png_file = '/Users/lding/Documents/MICA/outputs/GSE71585_Tasic/Tasic_k8_seurat_umap_true_label.png'
vs.scatter_plot(seurat_umap_true_label, out_png_file, method='umap', marker_size=2.0)









#%% PBMC
mica_dir = '/Users/lding/Documents/MICA/outputs'
ge_umap_out_file = '{}/PBMC_20k_default_parameters/40_1.4_misdist_0.6/clustering_umap_euclidean_40_1.4.txt'.format(mica_dir)
ge_tsne_out_file = '{}/PBMC_20k_default_parameters/40_1.4_misdist_0.6/clustering_tsne_euclidean_40_1.4.txt'.format(mica_dir)
mds_umap_out_file = '{}/PBMC_20k_MDS_new/cwl_lsf_k10_tsne_ClusterMem.txt'.format(mica_dir)
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                  'PBMC_sorted_20K_true_label.txt'
seurat_umap_out_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_umap.csv'


#%%
ge_umap = pd.read_csv(ge_umap_out_file, delimiter='\t', index_col=0)
ge_umap_true_label = copy.deepcopy(ge_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_umap_true_label['label'] = true_label

#%%
ge_tsne = pd.read_csv(ge_tsne_out_file, delimiter='\t', index_col=0)
ge_tsne_true_label = copy.deepcopy(ge_tsne)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_tsne_true_label['label'] = true_label

#%%
seurat_umap = pd.read_csv(seurat_umap_out_file, delimiter=',', index_col=0)
seurat_umap_true_label = copy.deepcopy(seurat_umap)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
seurat_umap_true_label['label'] = true_label
seurat_umap_true_label.columns = ['X', 'Y', 'label']

#%%
out_png_file = '{}/PBMC_20k_default_parameters/40_1.4_misdist_0.6/clustering_umap_euclidean_40_1.4_true_label.png'.format(mica_dir)
vs.scatter_plot(ge_umap_true_label, out_png_file, method='umap')

#%%
out_png_file = '{}/PBMC_20k_default_parameters/40_1.4_misdist_0.6/clustering_tsne_euclidean_40_1.4_true_label.png'.format(mica_dir)
vs.scatter_plot(ge_tsne_true_label, out_png_file, method='tsne')

#%%
out_png_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/seurat_umap_true_label.png'
vs.scatter_plot(seurat_umap_true_label, out_png_file, method='umap')

#%%
mds_umap = pd.read_csv(mds_umap_out_file, delimiter='\t', index_col=0)
mds_umap_true_label = copy.deepcopy(mds_umap)
mds_umap_true_label['label'] = true_label

#%%
out_png_file = '{}/PBMC_20k_MDS_new/cwl_lsf_k10_tsne_ClusterMem_true_label.png'.format(mica_dir)
vs.scatter_plot(mds_umap_true_label, out_png_file)






#%%
import seaborn as sns
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

#%%
summary = pd.read_excel('/Users/lding/Documents/MICA/outputs/summary.xlsx', index_col=0,
                        sheet_name='PBMC20k_100_neighbors')

#%%
summary = summary.pivot('dimension', 'resolution', 'ARI')

#%%
# sns.relplot(x="dimension", y="resolution", size="ARI", sizes=(15, 200), data=summary)
ax = plt.axes()
sns.heatmap(summary, cmap='coolwarm', ax=ax)
ax.set_title('Adjusted Rand index (PBMC 20k)')
plt.show()



#%%
summary = pd.read_excel('/Users/lding/Documents/MICA/outputs/summary.xlsx', index_col=0,
                        sheet_name='GSE60361_Ziesel')

#%%
summary = summary.pivot('dimension', 'resolution', 'ARI')

#%%
# sns.relplot(x="dimension", y="resolution", size="ARI", sizes=(15, 200), data=summary)
ax = plt.axes()
sns.heatmap(summary, cmap='coolwarm', ax=ax)
ax.set_title('Adjusted Rand index (Ziesel)')
plt.show()



#%%
summary = pd.read_excel('/Users/lding/Documents/MICA/outputs/summary.xlsx', index_col=0,
                        sheet_name='GSE71585_Tasic')

#%%
summary = summary.pivot('dimension', 'resolution', 'ARI')

#%%
# sns.relplot(x="dimension", y="resolution", size="ARI", sizes=(15, 200), data=summary)
ax = plt.axes()
sns.heatmap(summary, cmap='coolwarm', ax=ax)
ax.set_title('Adjusted Rand index (Tasic)')
plt.show()



#%%
summary = pd.read_excel('/Users/lding/Documents/MICA/outputs/summary.xlsx', index_col=0,
                        sheet_name='GSE75688_Chung')

#%%
summary = summary.pivot('dimension', 'resolution', 'ARI')

#%%
# sns.relplot(x="dimension", y="resolution", size="ARI", sizes=(15, 200), data=summary)
ax = plt.axes()
sns.heatmap(summary, cmap='coolwarm', ax=ax)
ax.set_title('Adjusted Rand index (Chung)')
plt.show()






#%%
summary = pd.read_excel('/Users/lding/Documents/MICA/Manuscript/Tables/ARI/ARI_GE.xlsx', index_col=0,
                        sheet_name='HMC_100_neighbors_all')

#%%
summary = summary.pivot('dimension', 'resolution', 'ARI')

#%%
# sns.relplot(x="dimension", y="resolution", size="ARI", sizes=(15, 200), data=summary)
ax = plt.axes()
sns.heatmap(summary, cmap='coolwarm', ax=ax)
ax.set_title('Adjusted Rand index (Human_Motor_Cortex)')
plt.show()








#%% PBMC CD4T
mica_dir = '/Users/lding/Documents/scMINER/PBMC20k_CD4'
ge_act_out_file = '{}/act_clustering_UMAP_euclidean_24_4.48169.txt'.format(mica_dir)
ge_exp_out_file = '{}/exp_clustering_UMAP_euclidean_24_3.32012.txt'.format(mica_dir)
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                  'PBMC_sorted_20K_true_label.txt'

#%%
ge_act = pd.read_csv(ge_act_out_file, delimiter='\t', index_col=0)
ge_act_true_label = copy.deepcopy(ge_act)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_act_true_label['label'] = true_label

#%%
out_png_file = '{}/act_clustering_UMAP_euclidean_24_4.48169.png'.format(mica_dir)
vs.scatter_plot(ge_act_true_label, out_png_file, method='umap')

#%%
ge_exp = pd.read_csv(ge_exp_out_file, delimiter='\t', index_col=0)
ge_exp_true_label = copy.deepcopy(ge_exp)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_exp_true_label['label'] = true_label

#%%
out_png_file = '{}/exp_clustering_UMAP_euclidean_24_3.32012.png'.format(mica_dir)
vs.scatter_plot(ge_exp_true_label, out_png_file, method='umap')

#%%
from sklearn.metrics.cluster import adjusted_rand_score
ge_act_true_label = copy.deepcopy(ge_act)
ge_act_true_label['true_label'] = true_label
ari_act = adjusted_rand_score(ge_act_true_label['label'], ge_act_true_label['true_label'])
print(ari_act)

#%%
ge_exp_true_label = copy.deepcopy(ge_exp)
ge_exp_true_label['true_label'] = true_label
ari_act = adjusted_rand_score(ge_exp_true_label['label'], ge_exp_true_label['true_label'])
print(ari_act)





#%% PBMC CD4T
mica_dir = '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_weightedmean_cell_type_specific_net/'
ge_act_out_file = '{}/clustering_UMAP_euclidean_24_8.16617.txt'.format(mica_dir)
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/SilverStd/PBMC_20k/' \
                  'PBMC_sorted_20K_true_label.txt'

#%%
ge_act = pd.read_csv(ge_act_out_file, delimiter='\t', index_col=0)
ge_act_true_label = copy.deepcopy(ge_act)
true_label = pd.read_csv(true_label_file, delimiter='\t', index_col=0)
ge_act_true_label['label'] = true_label

#%%
out_png_file = '{}/act_clustering_UMAP_euclidean_24_8.16617.png'.format(mica_dir)
vs.scatter_plot(ge_act_true_label, out_png_file, method='umap')

#%%
from sklearn.metrics.cluster import adjusted_rand_score
ge_act_true_label = copy.deepcopy(ge_act)
ge_act_true_label['true_label'] = true_label
ari_act = adjusted_rand_score(ge_act_true_label['label'], ge_act_true_label['true_label'])
print(ari_act)

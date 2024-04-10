#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#%%
sil_file = '/Users/lding/Documents/MICA/Manuscript/Tables/Silhouette/Silhouette_summary.xlsx'
sil_table = pd.read_excel(sil_file, index_col=0)
sil_table = sil_table.rename(columns={'Human_Motor_Cortex': 'H_Cortex'})

#%%
sub_table = sil_table.loc[['MICA', 'Seurat', 'SC3', 'Scanpy'],]
sub_table.loc['dataset'] = sub_table.columns

#%%
sub_table_T = sub_table.T

#%%
melt_sil_table = pd.melt(sub_table_T, id_vars=['dataset'])

#%%
melt_sil_table = melt_sil_table.rename(columns={'variable': 'method', 'value': 'Silhouette Index'})
melt_sil_table['Silhouette Index'] = melt_sil_table['Silhouette Index'].astype(float)

#%%
plt.close()
plt.figure(figsize=(6, 6))
sns.set(style="darkgrid", font='Arial')
ax = sns.barplot(x="method", y="Silhouette Index", hue="method", data=melt_sil_table, dodge=False, capsize=.2)
# plt.show()
plt.savefig('/Users/lding/Desktop/Silhouette_summary_pvalue.pdf', dpi=500)


#%%
import scipy
mica_sil = melt_sil_table.loc[melt_sil_table['method'] == 'MICA'].loc[:,'Silhouette Index']
seurat_sil = melt_sil_table.loc[melt_sil_table['method'] == 'Seurat'].loc[:,'Silhouette Index']
print(scipy.stats.wilcoxon(mica_sil, seurat_sil, alternative='greater'))


#%%
other_sil = melt_sil_table.loc[melt_sil_table['method'] != 'MICA'].loc[:,'Silhouette Index']
print(scipy.stats.ranksums(mica_sil, other_sil))
print(scipy.stats.mannwhitneyu(mica_sil, other_sil, alternative='greater'))

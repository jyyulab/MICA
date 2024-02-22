#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#%%
ARI_file = '/Users/lding/Documents/MICA/Manuscript/Tables/ARI/ARI_summary.xlsx'
ARI_table = pd.read_excel(ARI_file, index_col=0)
ARI_table = ARI_table.rename(columns={'Human_Motor_Cortex': 'H_Cortex'})

#%%
sub_table = ARI_table.loc[['MICA', 'Seurat', 'SC3', 'Scanpy'],]
sub_table.loc['dataset'] = sub_table.columns

#%%
sub_table_T = sub_table.T

#%%
melt_ARI_table = pd.melt(sub_table_T, id_vars=['dataset'])

#%%
melt_ARI_table = melt_ARI_table.rename(columns={'variable': 'method', 'value': 'ARI'})
melt_ARI_table['ARI'] = melt_ARI_table['ARI'].astype(float)

#%%
plt.close()
plt.figure(figsize=(6, 6))
sns.set(style="darkgrid", font='Arial')
ax = sns.barplot(x="method", y="ARI", hue="method", data=melt_ARI_table, dodge=False, capsize=.2)
# plt.show()
plt.savefig('/Users/lding/Desktop/ARI_summary_pvalue.pdf', dpi=500)


#%%
import scipy
mica_ari = melt_ARI_table.loc[melt_ARI_table['method'] == 'MICA'].loc[:,'ARI']
seurat_ari = melt_ARI_table.loc[melt_ARI_table['method'] == 'Seurat'].loc[:,'ARI']
print(scipy.stats.wilcoxon(mica_ari, seurat_ari, alternative='greater'))


#%%
other_ari = melt_ARI_table.loc[melt_ARI_table['method'] != 'MICA'].loc[:,'ARI']
print(scipy.stats.ranksums(mica_ari, other_ari))
print(scipy.stats.mannwhitneyu(mica_ari, other_ari, alternative='greater'))

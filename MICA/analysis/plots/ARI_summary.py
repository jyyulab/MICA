#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#%%
ARI_file = '/Users/lding/Documents/MICA/Manuscript/Tables/ARI/ARI_summary.xlsx'
ARI_table = pd.read_excel(ARI_file, index_col=0)

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
plt.figure(figsize=(10, 6))
sns.lineplot(data=melt_ARI_table, x="dataset", y="ARI", hue="method", style="method",
             markers=True, markersize=10)
# plt.show()

#%%
plt.savefig('/Users/lding/Desktop/ARI_summary.pdf', dpi=500)

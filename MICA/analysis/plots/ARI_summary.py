#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


#%%
ARI_file = '/Users/lding/Documents/MICA/Manuscript/Tables/ARI/ARI_summary.xlsx'
ARI_table = pd.read_excel(ARI_file, index_col=0, sheet_name='ARI_all_methods')
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
plt.figure(figsize=(12, 6))
sns.set(style="darkgrid", font='Arial')
# sns.lineplot(data=melt_ARI_table, x="dataset", y="ARI", hue="method", style="method",
#              markers=True, markersize=10)
# markers = {"MICA": "s", "Seurat": "s", "SC3": "s", "Scanpy": "s"}
sns.lineplot(data=melt_ARI_table, x="dataset", y="ARI", hue="method",
             marker='s')
# plt.show()

#%%
plt.savefig('/Users/lding/Desktop/AMI_summary_lineplot.pdf', dpi=500)






#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
# sns.set_theme(style="darkgrid", font='Arial')
sns.set(style="darkgrid", font='Arial')
g = sns.catplot(x='dataset', y='ARI', hue='method', data=melt_ARI_table, kind='bar', height=4, aspect=3.5)
plt.tight_layout()
plt.grid()
# plt.show()

#%%
plt.savefig('/Users/lding/Desktop/ARI_summary_catplot.pdf', dpi=500)





#%%
cell_count_table = ARI_table.loc[['MICA', 'Seurat', 'SC3', 'Scanpy'],]
cell_count_table.loc['#cells'] = ARI_table.loc['# cells']

#%%
cell_count_table_T = cell_count_table.T

#%%
melt_cell_count_table = pd.melt(cell_count_table_T, id_vars=['#cells'])

#%%
melt_cell_count_table = melt_cell_count_table.rename(columns={'variable': 'method', 'value': 'ARI'})
melt_cell_count_table['ARI'] = melt_cell_count_table['ARI'].astype(float)


#%%
plt.close()
# sns.set(style="ticks")
sns.set_theme(style="darkgrid", font='arial')
g2 = sns.catplot(x='#cells', y='ARI', hue='method', data=melt_cell_count_table, kind='bar', height=4, aspect=3.5)
plt.grid()

#%%
plt.savefig('/Users/lding/Desktop/ARI_summary_catplot_cell_count.pdf', dpi=500)

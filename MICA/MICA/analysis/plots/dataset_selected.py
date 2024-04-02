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
yan = melt_ARI_table.loc[melt_ARI_table['dataset'] == 'Yan']

#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
fig, ax = plt.subplots(figsize=(6, 8))
sns.set_theme(style="darkgrid", font='Arial')
sns.barplot(x='dataset', y='ARI', hue='method', data=yan, ax=ax)
plt.tight_layout()
# plt.grid()
# plt.show()

#%%
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .15)

#%%
plt.savefig('/Users/lding/Desktop/Yan_ARI_summary_barplot.pdf', dpi=500)






#%%
buettner = melt_ARI_table.loc[melt_ARI_table['dataset'] == 'Buettner']

#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
fig, ax = plt.subplots(figsize=(6, 8))
sns.set_theme(style="darkgrid", font='Arial')
sns.barplot(x='dataset', y='ARI', hue='method', data=buettner, ax=ax)
plt.tight_layout()
# plt.grid()
# plt.show()

#%%
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .15)

#%%
plt.savefig('/Users/lding/Desktop/Buettner_ARI_summary_barplot.pdf', dpi=500)






#%%
usoskin = melt_ARI_table.loc[melt_ARI_table['dataset'] == 'Usoskin']

#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
fig, ax = plt.subplots(figsize=(6, 8))
sns.set_theme(style="darkgrid", font='Arial')
sns.barplot(x='dataset', y='ARI', hue='method', data=usoskin, ax=ax)
plt.tight_layout()
# plt.grid()
# plt.show()

#%%
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .15)

#%%
plt.savefig('/Users/lding/Desktop/Usoskin_ARI_summary_barplot.pdf', dpi=500)






#%%
pbmc20k = melt_ARI_table.loc[melt_ARI_table['dataset'] == 'PBMC20k']

#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
fig, ax = plt.subplots(figsize=(6, 8))
sns.set_theme(style="darkgrid", font='Arial')
sns.barplot(x='dataset', y='ARI', hue='method', data=pbmc20k, ax=ax)
plt.tight_layout()
# plt.grid()
# plt.show()

#%%
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .15)

#%%
plt.savefig('/Users/lding/Desktop/PBMC20k_ARI_summary_barplot.pdf', dpi=500)




#%%
zeisel = melt_ARI_table.loc[melt_ARI_table['dataset'] == 'Zeisel']

#%%
plt.close()
# plt.figure(figsize=(15, 5))
# plt.rcParams["figure.figsize"] = (20, 5)
# plt.rcParams["xtick.labelsize"] = 7
# plt.rcParams["font.family"] = 'Arial'
# sns.set(rc={"figure.figsize": (20, 5)})
fig, ax = plt.subplots(figsize=(6, 8))
sns.set_theme(style="darkgrid", font='Arial')
sns.barplot(x='dataset', y='ARI', hue='method', data=zeisel, ax=ax)
plt.tight_layout()
# plt.grid()
# plt.show()

#%%
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

change_width(ax, .15)

#%%
plt.savefig('/Users/lding/Desktop/zeisel_ARI_summary_barplot.pdf', dpi=500)

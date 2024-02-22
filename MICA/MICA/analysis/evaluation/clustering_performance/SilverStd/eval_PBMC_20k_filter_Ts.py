#!/usr/bin/env python3
import pandas as pd


#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)


#%% Filter CD4+ T helper2 and CD8+ Cytotoxic T
filtered_cells = true_label[true_label['label'] != 'CD4+ T Helper2']

#%%
filtered_cells_again = filtered_cells[filtered_cells['label'] != 'CD8+ Cytotoxic T']

#%%
filtered_cells_again.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_16k.txt',
                            sep='\t')

#%%
from MICA.lib import preprocessing as pp
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'PBMC_20k'
input_file = '{}_MICA_input.h5ad'.format(author)
adata = pp.read_preprocessed_mat('{}/{}/{}/{}'.format(root_dir, level, author, input_file))


#%%
adata_filter = adata[filtered_cells_again['cell']]


#%%
adata_filter.write(filename='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/'
                            'PBMC_20k_MICA_input_filter_16k.h5ad')







#%% Filter CD4+ T helper2 and CD8+ Cytotoxic T
filtered_cells = true_label[true_label['label'] != 'CD4+ T Helper2']

#%%
filtered_cells_again = filtered_cells[filtered_cells['label'] != 'CD8+ Cytotoxic T']

#%%
filtered_cells_14k = filtered_cells_again[filtered_cells_again['label'] != 'CD34+']

#%%
filtered_cells_14k.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_14k.txt',
                          sep='\t')

#%%
adata_filter = adata[filtered_cells_14k['cell']]

#%%
adata_filter.write(filename='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/'
                            'PBMC_20k_MICA_input_filter_14k.h5ad')






#%% Filter CD4+ T helper2 and CD8+ Naive Cytotoxic T
filtered_cells = true_label[true_label['label'] != 'CD4+ T Helper2']

#%%
filtered_cells_again = filtered_cells[filtered_cells['label'] != 'CD8+/CD45RA+ Naive Cytotoxic']

#%%
filtered_cells_14k = filtered_cells_again[filtered_cells_again['label'] != 'CD34+']

#%%
filtered_cells_14k.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/filtered_cells_14k_2.txt',
                          sep='\t')

#%%
adata_filter = adata[filtered_cells_14k['cell']]

#%%
adata_filter.write(filename='/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/'
                            'PBMC_20k_MICA_input_filter_14k_2.h5ad')




#%%
root_dir = '/Users/lding/Documents/MICA/Datasets/HPC'
level = 'SilverStd'
author = 'PBMC_20k'
input_file = 'Filtered_DownSampled_SortedPBMC_data.csv'
pbmc20k_mat = pd.read_csv('{}/{}/{}/{}'.format(root_dir, level, author, input_file), header=0, index_col=0)


#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)


#%% Filter CD4+ T helper2 and CD8+ Cytotoxic T
filtered_cells = true_label[true_label['label'] != 'CD4+ T Helper2']

#%%
filtered_cells_again = filtered_cells[filtered_cells['label'] != 'CD8+ Cytotoxic T']

#%%
filtered_cells_14k = filtered_cells_again[filtered_cells_again['label'] != 'CD34+']

#%%
pbmc14k_mat = pbmc20k_mat.loc[filtered_cells_14k['cell'],:]

#%%
pbmc14k_mat.to_csv('/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/'
                   'PBMC_20k_MICA_input_filter_14k.txt', sep='\t')

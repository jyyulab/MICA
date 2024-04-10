#!/usr/bin/env python3

import pandas as pd
from collections import Counter


#%%
# mc_member_file = '/Users/lding/Documents/MICA/MetaCell/random_pivot/metacell_membership.txt'
# mc_member_file = '/Users/lding/Documents/MICA/MetaCell/closeness_pivot/metacell_membership.txt'
mc_member_file = '/Users/lding/Documents/MICA/MetaCell/closeness_pivot_top_20/metacell_membership.txt'
true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'

#%%
mc_member = pd.read_csv(mc_member_file, delimiter='\t', header=0)
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
mc_with_true_label = mc_member.merge(true_label, how='inner', left_on='barcode', right_on='cell')

#%%
impure_mc = []
num_mc = 0
mc_purity_total = 0.0
for l in set(mc_with_true_label.loc[:, 'clustering_label']):
    cells = mc_with_true_label.loc[mc_with_true_label['clustering_label'] == l]
    unique_labels = set(cells.loc[:, 'label'])
    if len(unique_labels) > 1:
        num_mc += 1
        impure_mc.append((l, unique_labels))
        cell_type_counter = Counter(cells.loc[:, 'label'])
        print(cell_type_counter)
        mc_purity_total += cell_type_counter.most_common(1)[0][1] / len(cells)

print(num_mc)
print(mc_purity_total / num_mc)

#%%
for mc in impure_mc:
    print(mc)
print(len(impure_mc))

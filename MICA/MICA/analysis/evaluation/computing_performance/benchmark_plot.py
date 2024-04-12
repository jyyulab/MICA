#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt


#%%
summary_df = pd.read_csv('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/summary.txt', sep='\t')

#%%
ax = sns.boxplot(x="num_neighbor_mi", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nnm_ari.pdf')

#%%
ax = sns.boxplot(x="num_neighbor_mi", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nnm_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="num_neighbor_mi", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nnm_vs_max_mem.pdf')





#%%
summary_df_nnm_100 = summary_df.loc[summary_df['num_neighbor_mi'] == 100]

#%%
ax = sns.boxplot(x="num_workers", y="run_time", data=summary_df_nnm_100, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nw_vs_run_time_nnm_100.pdf')

#%%
ax = sns.boxplot(x="num_workers", y="max_mem", data=summary_df_nnm_100, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nw_vs_max_mem_nnm_100.pdf')





#%%
ax = sns.boxplot(x="pruning_degree_multi", y="ARI", data=summary_df)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.15)

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/pdm_vs_ari_all.pdf')

#%%
ax = sns.boxplot(x="pruning_degree_multi", y="ARI", data=summary_df_nnm_100)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.15)

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/pdm_vs_ari_nnw_100.pdf')





#%%
summary_df = pd.read_csv('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/summary.txt', sep='\t')

#%%
ax = sns.boxplot(x="num_neighbor_euclidean", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/nne_vs_ari.pdf')

#%%
ax = sns.boxplot(x="num_neighbor_euclidean", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/nnm_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="num_neighbor_euclidean", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/nnm_vs_max_mem.pdf')





#%%
# summary_df_nnm_100 = summary_df.loc[summary_df['num_neighbor_mi'] == 100]

#%%
ax = sns.boxplot(x="num_workers", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/nw_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="num_workers", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne/nw_vs_max_mem.pdf')






#%%
summary_df = pd.read_csv('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/summary.txt', sep='\t')

#%%
ax = sns.boxplot(x="walk_length", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/walk_length_vs_ARI.pdf')

#%%
ax = sns.boxplot(x="num_walks", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/num_walks_vs_ARI.pdf')

#%%
ax = sns.boxplot(x="window_size", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/window_size_vs_ARI.pdf')



#%%
summary_df1 = summary_df.loc[summary_df['walk_length'] == 20]

#%%
summary_df2 = summary_df1.loc[summary_df1['window_size'] == 10]

#%%
ax = sns.boxplot(x="num_walks", y="ARI", data=summary_df2, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/num_walks_vs_ARI_wl_20_ws_10.pdf')




#%%
summary_df1 = summary_df.loc[summary_df['num_walks'] == 110]

#%%
summary_df2 = summary_df1.loc[summary_df1['window_size'] == 10]

#%%
ax = sns.boxplot(x="walk_length", y="ARI", data=summary_df2, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/walk_length_vs_ARI_nwalk_110_ws_10.pdf')




#%%
summary_df1 = summary_df.loc[summary_df['walk_length'] == 20]

#%%
summary_df2 = summary_df1.loc[summary_df1['num_walks'] == 110]

#%%
ax = sns.boxplot(x="window_size", y="ARI", data=summary_df2, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/window_size_vs_ARI_wl_20_nwalk_110.pdf')




#%%
summary_df = pd.read_csv('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/summary.txt', sep='\t')

#%%
ax = sns.boxplot(x="walk_length", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/walk_length_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="walk_length", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/walk_length_vs_max_mem.pdf')

#%%
ax = sns.boxplot(x="num_walks", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/num_walks_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="num_walks", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/num_walks_vs_max_mem.pdf')

#%%
ax = sns.boxplot(x="window_size", y="run_time", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/window_size_vs_run_time.pdf')

#%%
ax = sns.boxplot(x="window_size", y="max_mem", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/window_size_vs_max_mem.pdf')





#%%
summary_df1 = summary_df.loc[summary_df['walk_length'] == 20]

#%%
summary_df2 = summary_df1.loc[summary_df1['window_size'] == 10]

#%%
ax = sns.boxplot(x="num_walks", y="run_time", data=summary_df2, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/random_walk/num_walks_vs_run_time_wl_20_ws_10.pdf')





#%%
summary_df = pd.read_csv('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/hyperparameter/summary.txt', sep='\t')

#%%
ax = sns.boxplot(x="hyperparameter_p", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/hyperparameter/hyperp_vs_ARI.pdf')

#%%
ax = sns.boxplot(x="hyperparameter_q", y="ARI", data=summary_df, palette="coolwarm")

#%%
plt.savefig('/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/hyperparameter/hyperq_vs_ARI.pdf')


#%%
ax = sns.boxplot(x="hyperparameter_p", y="run_time", data=summary_df, palette="coolwarm")


#%%
ax = sns.boxplot(x="hyperparameter_q", y="run_time", data=summary_df, palette="coolwarm")


#%%
ax = sns.boxplot(x="hyperparameter_p", y="max_mem", data=summary_df, palette="coolwarm")


#%%
ax = sns.boxplot(x="hyperparameter_q", y="max_mem", data=summary_df, palette="coolwarm")

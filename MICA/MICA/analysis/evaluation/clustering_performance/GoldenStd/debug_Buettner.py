#!/usr/bin/env python3
import pandas as pd
import numpy as np
from sklearn import cluster
from MICA.bin import clustering
from sklearn.metrics.cluster import adjusted_rand_score
from multiprocessing import Pool
from functools import partial
from MICA.lib import utils
from MICA.lib import preprocessing
from MICA.lib import distance
import itertools
from sklearn import manifold


#%%
h5_file = '/private/tmp/docker_tmplg6f7ho7/Buettner_mi_0.h5'
hdf = pd.HDFStore(h5_file)
df1 = hdf['mi_0']

#%%
h5_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_dist.h5'
hdf = pd.HDFStore(h5_file)
print(hdf.keys())
df2 = hdf['mi']
df_norm = hdf['norm_mi']

#%%
h5_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_reduced.h5'
hdf = pd.HDFStore(h5_file)
print(hdf.keys())
df3 = hdf['mds']

#%%
def mds(in_mat_file, max_dim, out_file_name, perplexity=30, print_plot="True", dist_method="mi"):
    hdf = pd.HDFStore(in_mat_file)
    if dist_method == "mi":
        df = 1 - hdf["mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = df.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -H.dot(df ** 2).dot(H) / 2
    print(B)
    evals, evecs = utils.eigh(B, eigvals=(n - np.min([n, 200]), n - 1))
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    print(evals)
    evecs = evecs[:, idx]
    evals_pos = evals > 0
    print(evals_pos)
    L = np.diag(np.sqrt(evals[evals_pos]))
    print(L)
    V = evecs[:, evals_pos]
    Y = pd.DataFrame(
        data=V.dot(L),
        index=df.index,
        columns=["mds_" + str(x) for x in np.arange(1, L.shape[0] + 1)],
    )
    print(Y)

    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds

    if print_plot == "True":
        vis = utils.tsne(
            Y,
            max_dim,
            out_file_name,
            "mds",
            perplexity,
            print_plot,
        )
        vis.to_hdf(out_file_name + "_reduced", "mds_tsne")  # save preview in key "mds_tsne"


#%%
mds(h5_file, 19, 'Buettner_mds')

#%%
embedding = manifold.MDS(n_components=19, dissimilarity='precomputed')
X_transformed = embedding.fit_transform(df)


#%%
mi_file = '/Users/lding/Documents/MICA/tests/f8559b2d-3aad-43cc-a118-d727c72951cf/dimension_reduce/cwl_lsf_reduced.h5'
hdf = pd.HDFStore(mi_file)
df = hdf['mds']

#%%
n_cluster = 3
n_iter = 1
common_name = 'cwl_'

km = cluster.KMeans(n_clusters=3, max_iter=1000, n_init=1000)

km_results = pd.DataFrame(
    # data=np.transpose(km.fit_predict(df.iloc[:, 0:12])),
    data=np.transpose(km.fit_predict(X_transformed)),
    index=df.index,
    columns=["label"],
)


#%%
df_R = pd.read_csv("/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_R_MI.csv",
                   index_col=0, header=0)

#%%
km = cluster.KMeans(n_clusters=3, max_iter=1000, n_init=1000)

km_res = pd.DataFrame(
    data=np.transpose(km.fit_predict(df_R.iloc[:, 0:12])),
    index=df.index,
    columns=["label"],
)

#%%
mi_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_reduced.h5'
hdf = pd.HDFStore(mi_file)
df3 = hdf['mds']

km_results3 = pd.DataFrame(
    data=np.transpose(km.fit_predict(df3.iloc[:, 0:12])),
    index=df.index,
    columns=["label"],
)


#%%
def consensus_sc3(km_results, n_clusters, common_name=None):
    """ Implement SC3's consensus clustering. (https://www.nature.com/articles/nmeth.4236)
    Args:
        km_results (list of dataframes): each dataframe is a clustering results with two columns (cell_index,
                                         clustering_label)
        n_clusters (int): number of clusters
        common_name (name): common string to name the output files
    Returns:
        Clustering results
    """
    n_iter = len(km_results)
    if n_iter == 0:
        return None
    if len(km_results) == 0:
        return None

    conss_binary_mat = np.zeros((km_results[0].shape[0], km_results[0].shape[0]))
    for i in range(n_iter):
        arr = km_results[i].to_numpy()
        mask = arr[:, None] == arr
        binary_mat = mask[:,:,0].astype(int)
        conss_binary_mat += binary_mat
    conss_binary_mat = conss_binary_mat / n_iter

    clust = cluster.AgglomerativeClustering(
        linkage="ward", n_clusters=n_clusters, affinity="euclidean"
    )
    clust.fit(conss_binary_mat)

    cclust = pd.DataFrame(data=clust.labels_, index=km_results[0].index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})

    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    # conss_binary_mat.to_hdf(out_file_hdf, "membership")
    return cclust, out_file

#%%
def aggregate(km_results, n_clusters, common_name):
    n_iter = len(km_results)
    if n_iter == 0:
        return None
    mem = None
    for i in range(n_iter):
        df = km_results[i]
        dff = pd.DataFrame(data=df.values@df.T.values, index=df.index, columns=df.index)
        dff_div = pd.DataFrame(
            data=np.array((np.diag(dff),) * dff.shape[0]).T,
            index=dff.index,
            columns=dff.columns,
        )
        mem_mat = pd.DataFrame(
            data=dff / dff_div == 1,
            index=dff.index,
            columns=dff.columns,
            dtype=np.float32,
        )
        mem = mem_mat if i == 0 else mem + mem_mat.loc[mem.index, mem.columns]
        mem = mem / n_iter if i == n_iter - 1 else mem

    clust = cluster.AgglomerativeClustering(
        linkage="ward", n_clusters=n_clusters, affinity="euclidean"
    )
    clust.fit(mem)

    cclust = pd.DataFrame(data=clust.labels_, index=mem.index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})

    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    mem.to_hdf(out_file_hdf, "membership")
    return cclust, out_file


#%%
agg, out_f = consensus_sc3(km_results, n_cluster, common_name)

#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
# cluster_mem_file = '/Users/lding/Documents/MICA/outputs/GoldernStd/Buettner/mica-db8/cwl_lsf_k3_tsne_ClusterMem.txt'
# cluster_mem_file = '/Users/lding/Documents/MICA/outputs/GoldernStd/Buettner/' \
#                    'mica-f8559b2d-3aad-43cc-a118-d727c72951cf/cwl_lsf_k3_tsne_ClusterMem.txt'
cluster_mem_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/' \
                   'Buettner_k3_tsne_ClusterMem.txt'
# print(cluster_mem_file)
predict_label = pd.read_csv(cluster_mem_file, delimiter='\t', index_col=0)
# print(predict_label)
# agg.index = agg.index.astype(int)
# predict_label = agg
merged = true_label.merge(predict_label, left_on='cell', right_on='ID')
# print(merged)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print('{}'.format(ari))


#%%
true_label_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/yan/yan_true_label.txt'
true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)

#%%
# cluster_mem_file = '/Users/lding/Documents/MICA/outputs/GoldernStd/Buettner/mica-db8/cwl_lsf_k3_tsne_ClusterMem.txt'
# cluster_mem_file = '/Users/lding/Documents/MICA/outputs/GoldernStd/Buettner/' \
#                    'mica-f8559b2d-3aad-43cc-a118-d727c72951cf/cwl_lsf_k3_tsne_ClusterMem.txt'
cluster_mem_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/yan/' \
                   'Yan_k8_tsne_ClusterMem.txt'
# print(cluster_mem_file)
predict_label = pd.read_csv(cluster_mem_file, delimiter='\t', index_col=0)
predict_label.index = range(1, len(predict_label.index)+1)
predict_label.index.name = 'ID'
# print(predict_label)
# agg.index = agg.index.astype(int)
# predict_label = agg
merged = true_label.merge(predict_label, left_on='cell', right_on='ID')
# print(merged)
ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
print('{}'.format(ari))


#%%
adata = preprocessing.read_preprocessed_mat('/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/'
                                            'Buettner/Buettner_MICA_input.h5ad')
frame = adata.to_df()
ndarray = frame.to_numpy()

#%%
cell1 = ndarray[0,]
cell2 = ndarray[1,]
cell3 = ndarray[2,]

#%%
num_bins = int((ndarray.shape[0]) ** (1 / 3.0))
num_genes = ndarray.shape[1]

#%%
numba_dis = distance.numba_calc_mi_dis(cell1, cell2, num_bins, num_genes)
print(numba_dis)

#%%
print(distance.calc_norm_mi(cell1, cell2, num_bins, num_genes))
print(distance.calc_norm_mi(cell1, cell3, num_bins, num_genes))

#%%
print(utils.calc_mi(frame.iloc[0,], frame.iloc[1,], num_bins, num_genes))

#%%
norm_mi_mat = distance.calc_dis_mat(frame, frame, num_bins, num_genes)

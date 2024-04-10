#!/usr/bin/env python3
"""
This module contains helper functions, essential to the execution of MICA (Mutual Information-based
clustering algorithm).
"""

import sys
import umap as ump  # will pop up "joblib" deprecation warning message
import numpy as np
import pandas as pd
import matplotlib  # for plotting
matplotlib.use("Agg")
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from sklearn import cluster        # for kmeans
from sklearn import manifold       # for tsne
from sklearn import decomposition  # for PCA, etc
from scipy.cluster.hierarchy import dendrogram  # for heatmap
from scipy.linalg import eigh
from scipy.spatial import distance  # for euclidean distance


# def read_file(in_file, out_file_name):
#     """ Reads text file and stores data in a temporary HDF5-format file.
#
#     Args:
#         in_file_name  (str): path to input text file
#         out_file_name (str): user-defined name for output file
#     """
#     if in_file.endswith('.txt'):
#         frame = pd.read_csv(in_file, sep="\t", index_col=0).iloc[:, 0:]
#     if in_file.endswith('.h5ad'):
#         adata = anndata.read_h5ad(in_file)
#         frame = adata.to_df()
#     frame.to_hdf(out_file_name + ".h5.tmp", "slice_0")


def slice_file(df_file,  out_file_name, slice_size="1000"):
    """ Slices the HDF5 file.

    Determines the number of slices for the data, based on slice_size.
    Calculates start and end indices for slicing, creating a new dataframe
    based on those indices and appends them to the sliced output file using
    a unique-identifier in the format 'slice_00x' as key.

    Args:
        df_file       (str): path to HDF5-format file
        out_file_name (str): path to sliced output HDF5 file
        slice_size    (str): number of items in each slice
    """

    slice_size = int(slice_size)
    df = pd.HDFStore(df_file)["slice_0"]
    b = int(np.ceil(float(df.shape[0]) / float(slice_size)))
    digit = int(np.floor(np.log10(b)) + 1)
    for i in range(b):
        slice_name = str(i).zfill(digit)
        start = i * slice_size
        end = np.min([(i + 1) * slice_size, df.shape[0]])
        slice_ = pd.DataFrame(
            data=df.iloc[start:end, :].values,
            index=df.index[start:end],
            columns=df.columns,
        )
        slice_.to_hdf(out_file_name + ".sliced.h5", "slice_" + slice_name)
    pd.DataFrame(data=np.array(df.shape + (b,)), index=["row", "col", "slice"]).to_hdf(
        out_file_name + ".sliced.h5", "attr"
    )
    pd.DataFrame(data=df.columns, index=df.columns).to_hdf(
        out_file_name + ".sliced.h5", "cols"
    )
    pd.DataFrame(data=df.index, index=df.index).to_hdf(
        out_file_name + ".sliced.h5", "rows"
    )


def patch_file(df_file,  out_file_name):
    """ Prepares the HDF5 file for slicing. Completes the "temporary" HDF5-format file.

    Reads input file into several data frames. Indexes attributes as row, col and slice.
    Indexes columns and rows of data. Converts all data frames into an output HDF5 file
    with separate keys for each piece of data.
    
    Args:
        df_file       (str): path to HDF5-format file
        out_file_name (str): path to complete HDF5-format file
    """

    df = pd.HDFStore(df_file)["slice_0"]
    df.to_hdf(out_file_name + ".whole.h5", "slice_0")
    pd.DataFrame(data=np.array(df.shape + (1,)), index=["row", "col", "slice"]).to_hdf(
        out_file_name + ".whole.h5", "attr"
    )
    pd.DataFrame(data=df.columns, index=df.columns).to_hdf(
        out_file_name + ".whole.h5", "cols"
    )
    pd.DataFrame(data=df.index, index=df.index).to_hdf(
        out_file_name + ".whole.h5", "rows"
    )


def calc_prep(in_file, project_name):
    """ Prepares the already sliced input file for further calculation in MICA.
    
    Enters pairs of slices (matrices) into temporary HDF5-format files. It enters them
    individually, using their unique key. It also enters the parameter data for every single 
    pair into the key "params", which consists of: [key1, key2, num_bins, num_genes,
    pair_index, project_name, num_slices]
    
    Args:
        in_file      (str): path to sliced HDF5-format file
        project_name (str): project name used to generate path for final outputs
    """
    in_ = pd.HDFStore(in_file, "r")  # slice.h5
    # Use n^{1/3} as the bin size, where n is the number of genes.
    gene_count = int(in_["attr"].loc["col"])
    bins = int(np.floor(gene_count ** (1 / 3.0)))
    print('Number of genes: {}'.format(gene_count))
    print('Number of bins for estimating MI: {}'.format(bins))

    b = in_["attr"].loc["slice", 0]  # number of sliced matrix
    m = in_["attr"].loc["col", 0]  # number of genes
    digit = int(np.floor(np.log10(b)) + 1)  # some unique identifier
    total = int((b * (b + 1)) / 2)  # total number of calc jobs execute
    digit1 = int(np.floor(np.log10(total)) + 1)
    for i in range(b):
        key1 = "slice_" + str(i).zfill(digit)  # location of sliced matrix 1
        mat1 = in_[key1]  # sliced matrix 1def
        for j in range(i, b):
            key2 = "slice_" + str(j).zfill(digit)
            mat2 = in_[key2]
            idx = int(i * b + j - (i * (i + 1)) / 2)
            name = "mi_" + str(idx).zfill(digit1)  # name of MI pair
            mat_tmp = (project_name + "_" + name + ".h5.tmp")  # tmp directory as MIE_out/.tmp
            mat1.to_hdf(mat_tmp, key1)  # key: slice_00x
            mat2.to_hdf(mat_tmp, key2)
            pd.DataFrame(data=np.array([key1, key2, bins, m, name,  project_name, b]),
                         index=["key1", "key2",
                                "num_bins", "n_genes",
                                "MI_indx", "project_name", "num_slices"]).to_hdf(mat_tmp, "params")
    in_.close()


def vpearson(X, y):
    X_mean = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    y_mean = np.mean(y)
    r_num = np.sum((X-X_mean)*(y-y_mean), axis=1)
    r_den = np.sqrt(np.sum((X-X_mean)**2, axis=1)*np.sum((y-y_mean)**2))
    r = r_num/r_den
    return r


def calc_mi(arr1, arr2, bins, m):
    """ Calculates mutual information in between two cells, considering their gene expression levels
    
    This function is called by calc_distance_mat. It takes gene expression data from single cells,
    and compares them using standard calculation for mutual information. It builds a 2d histogram,
    which is used to calculate P(arr1, arr2)

    Args:
        arr1 (pandas series): gene expression data for a given cell in matrix_1
        arr2 (pandas series):
        bins           (int):
        m              (int):
    
    """
    fq = np.histogram2d(arr1.values, arr2.values, bins=(bins, bins))[0] / float(m)
    sm = np.sum(fq * float(m), axis=1)
    tm = np.sum(fq * float(m), axis=0)
    sm = np.asmatrix(sm / float(sm.sum()))
    tm = np.asmatrix(tm / float(tm.sum()))
    sm_tm = np.matmul(np.transpose(sm), tm)
    div = np.divide(fq, sm_tm, where=sm_tm != 0, out=np.zeros_like(fq))
    ent = np.log(div, where=div != 0, out=np.zeros_like(div))
    agg = np.multiply(fq, ent, out=np.zeros_like(fq), where=fq != 0)
    return agg.sum()


def calc_distance_mat(mat1, mat2, paras, method):
    """ Calculates a distance metric in between two matrices (slices)

    Calculates a distance metric using the preferred method of comparison. Iterates over each cell's
    gene expression data and populates a new matrix with rows and columns as cells from the input
    matrices. The resulting matrix is then converted to an HDF5-format file.

    Args:
        mat1  (pandas dataframe): a sliced part of the original matrix, with some fraction of the
                                  total cells as rows from original file and all gene expression
                                  attributes as columns
        mat2  (pandas dataframe): similar to mat1
        paras (pandas dataframe): a dataframe that holds an array of parameters from the whole dataset
        method             (str): the method to be used for the distance calculation (
                                        mutual information: "mi"
                                        euclidean distance: "euclidean"
                                        pearson correlation: "pearson"
                                        spearman correlation: "spearman")
    """

    bins = int(paras.loc["num_bins", 0])
    print('bins: {}'.format(bins))
    m = int(paras.loc["n_genes", 0])
    key = paras.loc["MI_indx", 0]

    project_name = paras.loc["project_name", 0]
    out_file_name = project_name + "_" + key + ".h5"
    print(out_file_name)

    if method == "mi":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index) 
        # start = time.time()
        for c in mat2.index:
            df.loc[mat1.index, c] = mat1.apply(
                calc_mi, axis=1, args=(mat2.loc[c, :], bins, m)
            )
        # end = time.time()
    elif method == "euclidean":
        dist = distance.cdist(mat1, mat2, method)
        df = pd.DataFrame(data=dist, index=mat1.index, columns=mat2.index, dtype="float")
    elif method == "pearson":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
        for c in mat2.index:
            df.loc[:, c] = vpearson(mat1.values, mat2.loc[c, :].values)
    elif method == "spearman":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
        for c in mat2.index:
            df.loc[:, c] = vpearson(mat1.rank(axis=1).values, mat2.loc[c, :].values)
    else:
        sys.exit("Distance Metrics not supported!\n")

    df.to_hdf(out_file_name, str(key))  # key: MI_indx
    paras.to_hdf(out_file_name, "params")


def merge_dist_mats(mi_slices, in_common_name, metrics):
    """ Iterates over and merges all distance matrices in one HDF5-format file

    Args:
        mi_slices    (str[]): a list of paths of sliced distance matrix files
        in_common_name (str): project name
        metrics        (str): metric used to calculate distance
    """
    mi_0 = pd.HDFStore(mi_slices[0])
    n_slice = int(mi_0["params"].loc["num_slices", 0])  # number of chopped dfs
    n = ((n_slice + 1) * n_slice) / 2

    if len(mi_slices) < n:
        sys.exit("Error - missing sliced MI file(s) for merging.")

    k = mi_0.keys()[0]
    mi_0[k].to_hdf(in_common_name + "_mi_whole.h5", key=k)
    mi_0.close()

    digit = int(np.floor(np.log10(n)) + 1)
    rows = []

    for i in range(len(mi_slices))[1:]:
        mi_k = pd.HDFStore(mi_slices[i])
        k = mi_k.keys()[0]
        mi_k[k].to_hdf(in_common_name + "_mi_whole.h5", key=k)
        mi_k.close()

    mi_whole = pd.HDFStore(in_common_name + "_mi_whole.h5")

    for i in range(n_slice):
        row_cols = []
        for j in range(i):
            idx = int(j * n_slice + i - (j * (j + 1)) / 2)
            key = "mi_" + str(idx).zfill(digit)  # mi_idx (name)
            mat_t = mi_whole[key]
            row_cols.append(mat_t.T)  # transposed MI add to merged file
        for j in range(i, n_slice):
            idx = int(i * n_slice + j - (i * (i + 1)) / 2)
            key = "mi_" + str(idx).zfill(digit)  # mi_idx (name)
            mat = mi_whole[key]
            row_cols.append(mat)
        rows.append(pd.concat(row_cols, axis=1))

    df = pd.concat(rows, axis=0)
    df.to_hdf(in_common_name + "_dist.h5", metrics)


def norm_mi_mat(in_mat_file, out_file_name):
    """Normalizes mutual information metric in the merged matrix
    
    Args:
        in_mat_file   (str): path to merged matrix
        out_file_name (str): name of output file
    """
    hdf = pd.HDFStore(in_mat_file)
    df = hdf["mi"]
    diag = np.asmatrix(np.diag(df))
    if in_mat_file == out_file_name + "_dist.h5":
        hdf.put("norm_mi", df / np.sqrt(np.matmul(diag.T, diag)))
    else:
        df.to_hdf(out_file_name + "_dist.h5", "mi")
        (df / np.sqrt(np.matmul(diag.T, diag))).to_hdf(
            out_file_name + "_dist.h5", "norm_mi"
        )
    hdf.close()


def tsne(data, max_dim, out_file_name, tag, perplexity=30, plot="True"):
    embed = manifold.TSNE(
        n_components=2,
        n_iter=5000,
        learning_rate=200,
        perplexity=perplexity,
        random_state=10,
        early_exaggeration=12.0).fit_transform(data.iloc[:, 0:max_dim])
    res = pd.DataFrame(data=embed, index=data.index, columns=["X", "Y"])
    if plot == "True":
        scatter(embed, out_file_name, tag)
    return res


def umap(data, max_dim, min_dist=0.25):
    embed = ump.UMAP(
        random_state=30,
        metric="euclidean",
        n_neighbors=10,
        min_dist=min_dist).fit_transform(data.iloc[:, 0:max_dim])
    res = pd.DataFrame(data=embed, index=data.index, columns=["X", "Y"])
    return res


def mds(in_mat_file, max_dim, out_file_name, perplexity=30, print_plot="True", dist_method="mi"):
    hdf = pd.HDFStore(in_mat_file)
    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = df.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -H.dot(df ** 2).dot(H) / 2
    evals, evecs = eigh(B, eigvals=(n - np.min([n, 200]), n - 1))
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    evals_pos = evals > 0
    L = np.diag(np.sqrt(evals[evals_pos]))
    V = evecs[:, evals_pos]
    Y = pd.DataFrame(
        data=V.dot(L),
        index=df.index,
        columns=["mds_" + str(x) for x in np.arange(1, L.shape[0] + 1)],
    )

    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds

    if print_plot == "True":
        vis = tsne(
            Y,
            max_dim,
            out_file_name,
            "mds",
            perplexity,
            print_plot,
        )
        vis.to_hdf(out_file_name + "_reduced", "mds_tsne")  # save preview in key "mds_tsne"


def lpl(in_mat_file, max_dim, out_file_name, perplexity=30, plot="True", dist_method="mi"):
    hdf = pd.HDFStore(in_mat_file)

    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = np.min(df.shape[0], 200)
    laplacian = manifold.SpectralEmbedding(
        n_components=n, eigen_solver="lobpcg", random_state=10
    )
    Y = pd.DataFrame(
        data=laplacian.fit_transform(df),
        index=df.index,
        columns=["lpl_" + str(x) for x in np.arange(1, n)],
    )
    Y.to_hdf(out_file_name + "_reduced.h5", "lpl")
    if plot == "True":
        tsne(Y, max_dim, out_file_name, "lpl", perplexity, plot)


def pca(in_mat_file, max_dim, out_file_name, perplexity=30, plot="True", dist_method="mi"):
    hdf = pd.HDFStore(in_mat_file)

    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = np.min(df.shape[0], 200)
    pca_ = decomposition.PCA(n_components=n, random_state=10)   # PCA on cell-cell distance matrix
                                                                # works for scRNA-seq data by SC3
    Y = pd.DataFrame(
        data=np.transpose(pca_.fit(df).components_),
        index=df.index,
        columns=["pca_" + str(x) for x in np.arange(1, n + 1)],
    )
    Y.to_hdf(out_file_name + "_reduced.h5", "pca")
    if plot == "True":
        tsne(Y, max_dim, out_file_name, "pca", perplexity, plot,)


def kmeans(in_mat, n_clusters, project_name, dim, bootstrap_id):
    out_file_name = project_name + "_kmeans_k" + str(n_clusters) + "_d" + str(dim) + ".h5.tmp." + str(bootstrap_id)
    km = cluster.KMeans(n_clusters=n_clusters, max_iter=1000, n_init=1000)

    km_res = pd.DataFrame(
        data=np.transpose(km.fit_predict(in_mat.iloc[:, 0:dim])),
        index=in_mat.index,
        columns=["label"],
    )
    # km_res.to_hdf(out_file_name, "kmeans")
    print("Executing kmeans:" + out_file_name)
    return km_res


# Deprecated, replaced by consensus.consensus_sc3
# def aggregate(km_results, n_clusters, common_name):
#     n_iter = len(km_results)
#     if n_iter == 0:
#         return None
#     mem = None
#     for i in range(n_iter):
#         df = km_results[i]
#         dff = pd.DataFrame(data=df.values@df.T.values, index=df.index, columns=df.index)
#         dff_div = pd.DataFrame(
#             data=np.array((np.diag(dff),) * dff.shape[0]).T,
#             index=dff.index,
#             columns=dff.columns,
#         )
#         mem_mat = pd.DataFrame(
#             data=dff / dff_div == 1,
#             index=dff.index,
#             columns=dff.columns,
#             dtype=np.float32,
#         )
#         mem = mem_mat if i == 0 else mem + mem_mat.loc[mem.index, mem.columns]
#         mem = mem / n_iter if i == n_iter - 1 else mem
#
#     clust = cluster.AgglomerativeClustering(
#         linkage="ward", n_clusters=n_clusters, affinity="euclidean"
#     )
#     clust.fit(mem)
#
#     cclust = pd.DataFrame(data=clust.labels_, index=mem.index, columns=["label"])
#     index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
#     label_map = {index[i]: i + 1000 for i in range(len(index))}
#     cclust = cclust.replace(to_replace={"label": label_map})
#     index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
#     map_back = {index[i]: i + 1 for i in range(len(index))}
#     cclust = cclust.replace(to_replace={"label": map_back})
#
#     out_file = common_name + "_k" + str(n_clusters)
#     out_file_hdf = out_file + "_cclust.h5"
#     cclust.to_hdf(out_file_hdf, "cclust")
#     mem.to_hdf(out_file_hdf, "membership")
#     return cclust, out_file


# def visualization(agg_mtx, reduced_mi_file, transformation, out_file_name, max_dim=0,
#                   visualize="umap", min_dist=0.25, perplexity=30,):
#     cclust = agg_mtx
#     hdf = pd.HDFStore(reduced_mi_file)
#     transformation = "pca" if transformation == "lpca" else transformation
#     hdf_trans = hdf[transformation.lower()]
#     hdf.close()
#
#     if visualize == "umap":
#         embed_2d = umap(hdf_trans, max_dim, min_dist)
#
#     elif visualize == "tsne":
#         perplexity = np.min([perplexity, np.max(cclust.groupby(["label"]).size())])
#         embed_2d = tsne(hdf_trans, max_dim, "", None, perplexity, "False")
#
#     cclust = pd.concat([cclust, embed_2d.loc[cclust.index, :]], axis=1)
#     res = cclust.loc[:, ["X", "Y", "label"]]
#     # save 2D embedding to txt file
#     out_file_name = out_file_name + "_" + visualize
#
#     scatter2(res, out_file_name + '.pdf')
#     res.to_csv(out_file_name + "_ClusterMem.txt", sep="\t")


def heatmap(df, linkage_, matrix, out_dir, out_file_name, tag, dendro="True"):
    fig = plt.figure(figsize=(16, 16), dpi=300)
    grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
    subgrid = gs.GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=grid[0],
        wspace=0.01,
        hspace=0.01,
        height_ratios=[1, 10],
        width_ratios=[10, 1],
    )
    if dendro == "True":
        ax = plt.Subplot(fig, subgrid[0])
        dendrogram(
            linkage_,
            ax=ax,
            orientation="top",
            labels=df.index,
            leaf_font_size=2,
            color_threshold=0,
        )
        ax.axis("off")
        fig.add_subplot(ax)
    ax1 = plt.Subplot(fig, subgrid[2])
    cax = ax1.matshow(matrix, cmap="YlOrRd", interpolation=None, aspect="auto")
    ax1.set_xticks([])
    ax1.set_yticks([])
    fig.add_subplot(ax1)
    ax3 = plt.Subplot(fig, subgrid[3])
    ax3.axis("off")
    cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
    plt.savefig(out_dir + out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


def heatmap2(data, labels, out_dir, out_file_name, tag):
    if data is None:
        return
    fig = plt.figure( figsize=(16, 16), dpi=300)
    grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
    subgrid = gs.GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=grid[0],
        wspace=0.001,
        hspace=0.1,
        height_ratios=[1, 20],
        width_ratios=[20, 1],
    )
    sorted_labels = labels.sort_values()
    bar = np.vstack((np.array(sorted_labels) + 1,) * 2) / np.max(np.unique(labels) + 1)
    ax1 = plt.Subplot(fig, subgrid[0])
    ax1.matshow(bar, cmap="jet", interpolation=None, aspect="auto", vmin=0, vmax=1)
    ax1.set_xticks([])
    ax1.set_yticks([])
    fig.add_subplot(ax1)
    ax2 = plt.Subplot(fig, subgrid[2])
    cax = ax2.matshow(
        data.loc[sorted_labels.index, sorted_labels.index],
        cmap="YlOrRd",
        interpolation=None,
        aspect="auto",
    )
    ax2.set_xticks([])
    ax2.set_yticks([])
    fig.add_subplot(ax2)
    ax3 = plt.Subplot(fig, subgrid[3])
    ax3.axis("off")
    cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
    plt.savefig(out_dir + out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


def scatter(df, out_file_name, tag, facecolor="none", edgecolor="r", marker="o", marker_size=20,):
    fig = plt.figure(figsize=(16, 16), dpi=300)
    plt.scatter(
        df[:, 0],
        df[:, 1],
        facecolor=facecolor,
        edgecolor=edgecolor,
        marker=marker,
        s=marker_size,
    )
    plt.ylabel("MICA-2")
    plt.xlabel("MICA-1")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


# def scatter2(data, out_file_name, marker_size=20, marker="o"):
#     if data is None:
#         return
#     fig = plt.figure(figsize=(10, 10), dpi=300)
#     lab = np.unique(data.loc[:, "label"])
#     colors = plt.cm.jet(np.linspace(0, 1, len(lab)))
#     for z in lab:
#         df = data.loc[data.loc[:, "label"] == z, ["X", "Y"]]
#         plt.scatter(
#             df.loc[:, "X"],
#             df.loc[:, "Y"],
#             facecolor=colors[z-1],
#             s=marker_size,
#             marker=marker,
#             vmin=0,
#             vmax=len(lab),
#             label=str(z) + "(" + str(df.shape[0]) + ")",
#             alpha=0.7,
#         )
#         center = np.mean(df, axis=0)
#         plt.scatter(
#             center.loc["X"],
#             center.loc["Y"],
#             marker="o",
#             c="white",
#             alpha=0.7,
#             s=100,
#             edgecolor="k",
#         )
#         plt.scatter(
#             center.loc["X"],
#             center.loc["Y"],
#             marker="$%d$" % z,
#             c="black",
#             alpha=0.7,
#             s=80,
#             edgecolor="k",
#         )
#     plt.ylabel("MICA-2")
#     plt.xlabel("MICA-1")
#     plt.xticks([])
#     plt.yticks([])
#     plt.legend(
#         loc="center left",
#         bbox_to_anchor=(1, 0.5),
#         title="Clusters(" + str(data.shape[0]) + ")",
#     )
#     plt.savefig(
#         out_file_name,
#         bbox_inches="tight",
#     )

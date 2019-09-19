#!/usr/bin/env python3
import sys
import itertools
import time
from multiprocessing import Pool
from functools import partial
import pandas as pd
import utils

start_time = time.time()


def km_multiprocess(
    mi_file,
    n_cluster,
    n_iter,
    common_name,
    dims=[19],
    thread=1,
):
    pool = Pool(processes=thread)
    hdf = pd.HDFStore(mi_file)
    r = []

    for trans in hdf.keys():
        df = hdf[trans]

        # def utils.kmeans(in_mat, n_clusters, project_name, dim, bootstrap_id)
        km_iterable = partial(utils.kmeans, df, n_cluster, common_name)
        iterations = itertools.product(dims, range(n_iter))

        res = pool.starmap(km_iterable, iterations)
        r = r + res

    pool.close()
    hdf.close()
    return r


if __name__ == '__main__':

    in_file = sys.argv[1]
    dr = sys.argv[2].lower()
    k = int(sys.argv[3])
    n_bootstrap = int(sys.argv[4])
    out_name = sys.argv[5]

    plot_method = sys.argv[6]
    umap_min_dist = float(sys.argv[7])
    tsne_perplexity = sys.argv[8]
    plot_dim = int(sys.argv[9])

    n_thread = int(sys.argv[10])
    dim_km = sys.argv[11:]
    dim_km = map(int, dim_km)

    # in_file = "/home/cqian/PBMC12K/PBMC12k_reduced.h5"
    # dr = "mds"
    # k = 8
    # n_bootstrap = 10
    # out_name = "pbmc_psub"
    # dim_km = [19]
    # n_thread = 6
    # plot_method = "tsne"
    # plot_dim = 19
    # umap_min_dist = 0.1
    # tsne_perplexity = 30

    result = km_multiprocess(in_file, n_cluster=k, n_iter=n_bootstrap,
                             common_name=out_name, dims=dim_km, thread=n_thread)

    # def aggregate(result, n_clusters, common_name):
    agg, out_f = utils.aggregate(result, k, out_name)

    utils.visualization(agg,  # consensus clustering result
                        in_file,  # reduced_mi_file
                        dr,  # transformation
                        out_f,  # output file name
                        max_dim=plot_dim,
                        visualize=plot_method.lower(),
                        min_dist=float(umap_min_dist),
                        perplexity=int(tsne_perplexity)
                        )

    print("--- %s seconds ---" % (time.time() - start_time))
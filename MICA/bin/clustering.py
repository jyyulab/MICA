#!/usr/bin/env python3

import itertools
import time
import argparse
import pandas as pd
from multiprocessing import Pool
from functools import partial
from MICA.lib import utils


def main():
    """Handles arguments and calls the driver function"""
    head_description = "Clusters reduced data using consensus k-mean"
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Reduced matrix file in h5 format')
    parser.add_argument('-dr', '--merge-option', metavar='STR', required=True,
                        help='Method used for dimension reduction')
    parser.add_argument('-k', '--k', type=int, metavar='INT', required=True,
                        help='k value for k-means')
    parser.add_argument('-n', '--n-bootstrap', type=int, metavar='INT', required=True,
                        help='Number of bootstraps')
    parser.add_argument('-o', '--output-file', metavar='STR', required=True,
                        help='Output file name')
    parser.add_argument('-m', '--plot-method', metavar='STR', required=True,
                        help='Method for plotting after clustering')
    parser.add_argument('-d', '--umap-min-dist', type=float, metavar='FLOAT', required=True,
                        help='Minimum distance for umap')
    parser.add_argument('-p', '--tsne-perplexity', type=int, metavar='INT', required=True,
                        help='TSNE perplexity measure')
    parser.add_argument('-dim', '--plot-dim', type=int, metavar='INT', required=True,
                        help='Dimension of plot')
    parser.add_argument('-t', '--n-thread', type=int, metavar='INT', required=True,
                        help='Total number of threads used')
    parser.add_argument('-km', '--dim-km', nargs='+', type=int, metavar='INT', required=True,
                        help='Size of k-means dimensions')
    args = parser.parse_args()
    
    clustering(args.input_file, args.merge_option.lower(), args.k, args.n_bootstrap, args.output_file, args.plot_method,
               args.umap_min_dist, args.tsne_perplexity, args.plot_dim, args.n_thread, args.dim_km)


def km_multiprocess(mi_file, n_cluster, n_iter, common_name, dims=[19], num_processes=1):
    pool = Pool(processes=num_processes)
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


def clustering(in_file, dr, k, n_bootstrap, out_name,
               plot_method, umap_min_dist, tsne_perplexity, plot_dim, n_processes, dim_km):
    start_time = time.time()
    dim_km = map(int, dim_km)
    result = km_multiprocess(in_file, n_cluster=k, n_iter=n_bootstrap,
                             common_name=out_name, dims=dim_km, num_processes=n_processes)

    # def aggregate(result, n_clusters, common_name):
    agg, out_f = utils.aggregate(result, k, out_name)

    utils.visualization(agg,  # consensus clustering result
                        in_file,  # reduced_mi_file
                        dr,  # transformation
                        out_f,  # output file name
                        max_dim=plot_dim,
                        visualize=plot_method.lower(),
                        min_dist=umap_min_dist,
                        perplexity=tsne_perplexity
                        )
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    main()

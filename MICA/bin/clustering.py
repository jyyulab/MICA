#!/usr/bin/env python3

import logging
import itertools
import argparse
import pandas as pd
from multiprocessing import Pool
from functools import partial
from MICA.lib import utils
from MICA.lib import consensus
from MICA.lib import clustering as cl
from MICA.lib import visualize as vi


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

        km_iterable = partial(utils.kmeans, df, n_cluster, common_name)
        iterations = itertools.product(dims, range(n_iter))

        res = pool.starmap(km_iterable, iterations)
        r = r + res

    pool.close()
    hdf.close()
    return r


def clustering(in_file, dr, k, n_bootstrap, out_name,
               plot_method, umap_min_dist, tsne_perplexity, plot_dim, n_processes, dim_km):
    dim_km = map(int, dim_km)
    result = km_multiprocess(in_file, n_cluster=k, n_iter=n_bootstrap,
                             common_name=out_name, dims=dim_km, num_processes=n_processes)

    agg, out_f = consensus.consensus_sc3(result, k, out_name)

    hdf = pd.HDFStore(in_file)
    dr = "pca" if dr == "lpca" else dr
    hdf_trans = hdf[dr.lower()]
    hdf.close()

    # UMAP or tSNE plot
    embed_res = vi.visual_embed(agg, hdf_trans, '.', dr_dim=19, visual_method=plot_method.lower(),
                                min_dist=umap_min_dist, perplexity=tsne_perplexity, marker_size=5.0,
                                marker_scale=5.0)

    # UMAP or tSNE plot
    # utils.visualization(agg,        # consensus clustering result
    #                     in_file,    # reduced_mi_file
    #                     dr,         # transformation
    #                     out_f,      # output file name
    #                     max_dim=plot_dim,
    #                     visualize=plot_method.lower(),
    #                     min_dist=umap_min_dist,
    #                     perplexity=tsne_perplexity)
    #
    # Silhouette plot
    # hdf = pd.HDFStore(in_file)
    # dr = "pca" if dr == "lpca" else dr
    # hdf_trans = hdf[dr.lower()]
    # hdf.close()

    silhouette_avg = cl.silhouette_analysis(embed_res, k, embed_res.loc[:, ['X', 'Y']], '.')
    logging.info('Silhouette score: {}'.format(silhouette_avg))


if __name__ == '__main__':
    main()

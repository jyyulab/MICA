#!/usr/bin/env python3
# MICA Version that uses a graph embedding method for dimension reduction on MI-kNN graph.

import sys
import logging
import time
import argparse
import numpy as np
import pandas as pd
import pathlib
from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs
from MICA.lib import consensus as cs


def main():
    head_description = 'MICA - Mutual Information-based Clustering Analysis tool. This version uses a graph ' \
                       'embedding method for dimension reduction on MI-kNN graph.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                        help='Path to an input file (h5ad file or tab-delimited text file)')
    parser.add_argument('-o', '--output-dir', metavar='DIR', required=True, help='Path to final output directory')
    parser = add_ge_arguments(parser)
    parser.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                        help='Visualization method UMAP or t-SNE (default: UMAP)')
    parser.add_argument('-md', '--min-dist', metavar='FLOAT', required=False, default=0.6, type=float,
                        help='min_dist parameter in UMAP, minimum distance of points in the embedded space '
                             '(default: 0.6)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mica_ge(args)


def add_ge_arguments(parser):
    parser.add_argument('-dm', '--dr-method', metavar='STR', required=False, choices=['node2vec', 'deepwalk'],
                        default='node2vec', help='Dimension reduction method [node2vec | deepwalk] (default: node2vec)')
    parser.add_argument('-dd', '--dr-dim', metavar='INT', required=False, default=20, type=int,
                        help='Number of dimensions to reduce to (default: 20)')
    parser.add_argument('-ir', '--min-resolution', metavar='FLOAT', required=False, default=0.2, type=float,
                        help='Determines the minimum size of the communities (default: 0.2)')
    parser.add_argument('-ar', '--max-resolution', metavar='FLOAT', required=False, default=4.2, type=float,
                        help='Determines the maximum size of the communities (default: 4.2)')
    parser.add_argument('-ss', '--step-size', metavar='FLOAT', required=False, default=0.4, type=float,
                        help='Determines the step size to sweep resolution from min_resolution to max_resolution '
                             '(default: 0.4)')
    parser.add_argument('-nw', '--num-workers', metavar='INT', required=False, default=10, type=int,
                        help='Number of works to run in parallel (default: 10)')
    return parser


def mica_ge(args):
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    adata = pp.read_preprocessed_mat(args.input_file)
    frame = adata.to_df()
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building MI-based kNN graph ...')
    knn_indices, knn_dists = ng.nearest_neighbors_NNDescent(frame.to_numpy(), num_jobs=args.num_workers)
    knn_graph = ng.build_graph_from_indices(knn_indices, knn_dists)

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    edgelist_file = '{}/NNDescent_knn_graph_edgelist.txt'.format(args.output_dir)
    with open(edgelist_file, 'w') as fout:
        for edge in knn_graph.edges():
            fout.write('{}\t{}\t{}\n'.format(edge[0], edge[1], knn_graph.get_edge_data(edge[0], edge[1])['MI']))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Performing dimension reduction using {} method ...'.format(args.dr_method))
    emb_file = '{}/knn_graph_emb_{}_{}.txt'.format(args.output_dir, args.dr_method, args.dr_dim)
    if args.dr_method == 'node2vec':
        dr.dim_reduce_node2vec_pecanpy(edgelist_file, emb_file, dim=args.dr_dim, num_jobs=args.num_workers)
        # wv = dr.dim_reduce_node2vec(knn_graph, dim=args.dr_dim, walk_len=10, n_walks=10)
        # print(wv)
    elif args.dr_method == 'deepwalk':
        # dr.dim_reduce_deepwalk(edgelist_file, emb_file, dim=args.dr_dim)
        sys.exit('Error - deepwalk has not been tested: {}'.format(args.dr_method))
    else:
        sys.exit('Error - invalid dimension reduction method: {}'.format(args.dr_method))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Performing clustering ...')
    mat_dr_df = pd.read_csv(emb_file, delimiter=' ', skiprows=1, index_col=0, names=np.arange(1, args.dr_dim+1))
    mat_dr_df.sort_index(inplace=True)
    mat_dr = mat_dr_df.to_numpy()
    logging.info(mat_dr.shape)
    G = ng.build_graph(mat_dr, dis_metric='euclidean')
    partitions = cl.graph_clustering_parallel(G, min_resolution=args.min_resolution, max_resolution=args.max_resolution,
                                              step_size=args.step_size, num_workers=args.num_workers)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Group partitions based on number of clusters and consensus clustering ...')
    cluster_dict = cs.group_partition(partitions, frame.index)
    aggs = []
    for num_cluster in cluster_dict.keys():
        agg, out_f = cs.consensus_sc3(cluster_dict[num_cluster], num_cluster, 'consensus')
        aggs.append((agg, num_cluster))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('Visualizing clustering results using {}'.format(args.visual_method))
    for agg, num_cluster in aggs:
        logging.info('Number of clusters: {}'.format(num_cluster))
        vs.visual_embed(agg, mat_dr, args.output_dir, suffix=num_cluster,
                        visual_method=args.visual_method, num_works=args.num_workers, min_dist=args.min_dist)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    out_h5_file = '{}/clustered.h5ad'.format(args.output_dir)
    logging.info('Write h5ad file {} ...'.format(out_h5_file))
    pp.write_h5(adata, out_h5_file)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


if __name__ == "__main__":
    main()

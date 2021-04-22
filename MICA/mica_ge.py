#!/usr/bin/env python3
# MICA Version that uses a graph embedding method for dimension reduction on MI-kNN graph.

import sys
import logging
import time
import argparse
import numpy as np
import pandas as pd
from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs


def main():
    head_description = 'MICA is a Mutual Information-based nonlinear Clustering Analysis tool designed for \
scRNA-seq data. This version uses a graph embedding method for dimension reduction on MI-kNN graph.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                        help='Path to an input file (h5ad file or tab-delimited text file)')
    parser.add_argument('-o', '--output-dir', metavar='DIR', required=True, help='Path to final output directory')
    parser.add_argument('-m', '--dr-method', metavar='STR', required=False, choices=['node2vec', 'deepwalk'],
                        default='node2vec', help='Dimension reduction method [node2vec | deepwalk] (default: node2vec)')
    parser.add_argument('-d', '--dr-dim', metavar='INT', required=False, default=12, type=int,
                        help='Number of dimensions to reduce to (default: 20)')
    parser.add_argument('-e', '--resolution', metavar='FLOAT', required=False, default=1.0, type=float,
                        help='Determines size of the communities. (default: 1.4)')
    parser.add_argument('-v', '--visual-method', metavar='STR', required=False, default='umap', type=str,
                        choices=['umap', 'tsne'], help='Visualization embedding method [umap | tsne] (default: umap)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    adata = pp.read_preprocessed_mat(args.input_file)
    frame = adata.to_df()
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building MI-based kNN graph ...')
    knn_indices, knn_dists = ng.nearest_neighbors_NNDescent(frame.to_numpy())
    knn_graph = ng.build_graph_from_indices(knn_indices, knn_dists)
    edgelist_file = '{}/knn_graph.edgelist.txt'.format(args.output_dir)
    with open(edgelist_file, 'w') as fout:
        for edge in knn_graph.edges():
            fout.write('{}\t{}\t{}\n'.format(edge[0], edge[1], knn_graph.get_edge_data(edge[0], edge[1])['MI']))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Performing dimension reduction using {} method ...'.format(args.dr_method))
    emb_file = '{}/knn_graph.{}.emb.txt'.format(args.output_dir, args.dr_method)
    if args.dr_method == 'node2vec':
        dr.dim_reduce_node2vec_pecanpy(edgelist_file, emb_file, dim=args.dr_dim)
        # wv = dr.dim_reduce_node2vec(knn_graph, dim=args.dr_dim, walk_len=10, n_walks=10)
        # print(wv)
    elif args.dr_method == 'deepwalk':
        dr.dim_reduce_deepwalk(edgelist_file, emb_file, dim=args.dr_dim)
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
    G = ng.build_graph(mat_dr, dis_metric='euclidean')
    partition = cl.graph_clustering(G, resolution=args.resolution)
    # partition = cl.graph_clustering(knn_graph, resolution=args.resolution)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Visualizing clustering results using {} ...'.format(args.visual_method))
    vs.visual_embed_louvain(partition, args.resolution, frame.index, mat_dr, args.output_dir,
                            embed_method=args.visual_method)
    # vs.visual_embed_louvain(partition, args.resolution, frame.index, frame, args.output_dir,
    #                         embed_method=args.visual_method)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    out_h5_file = '{}/clustered.h5ad'.format(args.output_dir)
    logging.info('Write h5ad file {} ...'.format(out_h5_file))
    pp.write_h5(adata, out_h5_file, partition, args.resolution)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


if __name__ == "__main__":
    main()

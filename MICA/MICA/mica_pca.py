#!/usr/bin/env python3

import sys
import argparse
import time
import logging
import scipy
import numpy as np
from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs


def main():
    head_description = '''MICA is a mutual information-based nonlinear clustering analysis tool designed for 
    scRNA-seq data. This version uses PCA for dimension reduction on euclidean space.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)
    parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                        help='Path to an input file (h5ad file or tab-delimited text file)')
    parser.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                        help='Path to final output directory')
    parser.add_argument('-dm', '--dr-method', metavar='STR', required=False, default='pca', type=str,
                        help='Dimension reduction method [pca | mds ] (default: pca)')
    parser.add_argument('-dd', '--dr-dim', metavar='INT', required=False, default=50, type=int,
                        help='Number of dimensions to reduce features of the input matrix to')
    parser.add_argument('-nn', '--num-neighbors', metavar='INT', required=False, default=20, type=int,
                        help='Number of neighbors of building neighboring graph (default: 20)')
    parser.add_argument('-ir', '--min-resolution', metavar='FLOAT', required=False, default=0.2, type=float,
                        help='Determines the minimum size of the communities (default: 0.2)')
    parser.add_argument('-ar', '--max-resolution', metavar='FLOAT', required=False, default=3.4, type=float,
                        help='Determines size of the communities. (default: 3.4)')
    parser.add_argument('-ss', '--step-size', metavar='FLOAT', required=False, default=0.4, type=float,
                        help='Determines the step size to sweep resolution from min_resolution to max_resolution '
                             '(default: 0.4)')
    parser.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                        choices=['UMAP', 't-SNE'], help='Visualization embedding method [UMAP | t-SNE] (default: UMAP)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    adata = pp.read_preprocessed_mat(args.input_file)
    print(adata.shape)
    if type(adata.X) is scipy.sparse.csr.csr_matrix:
        ndarr = adata.X.toarray()
    elif type(adata.X) is np.ndarray:
        ndarr = adata.X
    else:
        sys.exit('Error - invalid adata.X type: {}'.format(type(adata.X)))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Dimension reduction to {} dimensions using {} method ...'.format(args.dr_dim, args.dr_method))
    frame_dr = dr.dim_reduce_global(ndarr, dim=args.dr_dim)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building {}-nearest neighbor graph ...'.format(args.num_neighbors))
    G = ng.build_graph(frame_dr, num_neighbors=args.num_neighbors)
    # nx.write_gml(G, '{}/knn_{}_{}.gml'.format(args.output_dir, args.dist_metric, args.num_neighbors))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Performing graph clustering ...')
    partitions = cl.graph_clustering(G, max_resolution=args.max_resolution)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Visualizing clustering results using {} ...'.format(args.visual_method))
    # for partition in partitions:
    #     vs.visual_embed(partition, adata.obs.index, frame_dr, args.output_dir, visual_method=args.visual_method)

    for i, resolution in enumerate(list(np.arange(args.min_resolution, args.max_resolution+0.1, args.step_size))):
        resolution_round = np.round(resolution, 2)
        logging.info('Louvain resolution: {}'.format(resolution_round))
        vs.visual_embed(partitions[i], adata.obs.index, frame_dr, args.output_dir, resolution=resolution_round,
                        visual_method=args.visual_method)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


if __name__ == "__main__":
    main()

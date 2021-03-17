#!/usr/bin/env python3

import sys
import argparse
import time
import logging
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
    parser.add_argument('-r', '--dr-method', metavar='STR', required=False, default='pca', type=str,
                        help='Dimension reduction method [pca | mds ] (default: pca)')
    parser.add_argument('-d', '--dr-dim', metavar='INT', required=False, default=50, type=int,
                        help='Number of dimensions to reduce features of the input matrix to')
    parser.add_argument('-n', '--num-neighbors', metavar='INT', required=False, default=20, type=int,
                        help='Number of neighbors of building neighboring graph (default: 20)')
    parser.add_argument('-e', '--resolution', metavar='FLOAT', required=False, default=1.0, type=float,
                        help='Determines size of the communities. (default: 1.0)')
    parser.add_argument('-v', '--visual-method', metavar='STR', required=False, default='umap', type=str,
                        choices=['umap', 'tsne'], help='Visualization embedding method [umap | tsne] (default: umap)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    frame = pp.read_preprocessed_mat(args.input_file)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Dimension reduction to {} dimensions using {} method ...'.format(args.dr_dim, args.dr_method))
    frame_dr = dr.dim_reduce_global(frame, dim=args.dr_dim)
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
    partition = cl.graph_clustering(G, resolution=args.resolution)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Visualizing clustering results using {} ...'.format(args.visual_method))
    vs.visual_embed(partition, frame.index, frame_dr, args.output_dir, embed_method=args.visual_method)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


if __name__ == "__main__":
    main()

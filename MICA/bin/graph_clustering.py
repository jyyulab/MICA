#!/usr/bin/env python3

import argparse
import pandas as pd
from sklearn.neighbors import kneighbors_graph
import networkx as nx
import community
from MICA.lib import utils


def main():
    """Handles arguments and calls the driver function."""
    head_description = "Graph-based clustering of reduced matrix."
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Reduced matrix file in h5 format')
    parser.add_argument('-dr', '--dim-reduce', metavar='STR', required=True,
                        help='Method used for dimension reduction')
    parser.add_argument('-k', '--k', type=int, metavar='INT', required=True,
                        help='k in k nearest neighbor')
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
    parser.add_argument('-ds', '--dims', type=int, metavar='INT', required=True,
                        help='Number of reduced dimensions to be used for clustering')
    args = parser.parse_args()

    G = build_kNN(args.input_file, args.dim_reduce, args.k, args.dims)
    louvein_clustering(G, args.input_file, args.dim_reduce, args.output_file,
                       args.plot_method, args.umap_min_dist, args.tsne_perplexity, args.plot_dim)


def build_kNN(in_file, dr, num_neighbors, dim=19):
    """ Build k nearest neighbour graph of cells based on a dimension reduced distance matrix.
    Args:
        in_file (file): Reduced matrix file in h5 format
        dr (STR): Method used for dimension reduction
        num_neighbours (STR): number of neighbours
        dim (int): number of dimensions to be used for the clustering
    """
    hdf = pd.HDFStore(in_file)
    df = hdf[dr]
    kneighbor_graph = kneighbors_graph(df.iloc[:, 0:dim], n_neighbors=num_neighbors, include_self=False).toarray()
    G = nx.from_numpy_matrix(kneighbor_graph)
    return G


def louvein_clustering(G, in_file, dr, out_name, plot_method, umap_min_dist, tsne_perplexity, plot_dim):
    """ Perform louvien clustering on kNN, then visualize. """
    partition = community.best_partition(G)
    hdf = pd.HDFStore(in_file)
    index = hdf[dr].index
    labels = [x + 1 for x in partition.values()]
    cclust = pd.DataFrame(data=labels, index=index, columns=["label"])

    utils.visualization(cclust,  # consensus clustering result
                        in_file,  # reduced_mi_file
                        dr,  # transformation
                        out_name,  # output file name
                        max_dim=plot_dim,
                        visualize=plot_method.lower(),
                        min_dist=umap_min_dist,
                        perplexity=tsne_perplexity
                        )


if __name__ == '__main__':
    main()

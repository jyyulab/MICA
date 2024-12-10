#!/usr/bin/env python3
# MICA GE version uses a graph embedding method for dimension reduction on MI-kNN graph.

import sys
import logging
import time
import argparse
import numpy as np
import pandas as pd
import pathlib
import networkx as nx
import hnswlib
from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs
from MICA.lib import consensus as cs
from MICA.lib import metacell as mc


np.set_default_dtype(np.float32)


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
    parser.add_argument('-da', '--dr-modality', metavar='STR', required=False, choices=['gene', 'cell'],
                        default='gene', help='Dimension reduction modality [gene | cell] (default: gene)')
    parser.add_argument('-dm', '--dr-method', metavar='STR', required=False, choices=['node2vec', 'deepwalk'],
                        default='node2vec', help='Dimension reduction method [node2vec | deepwalk] (default: node2vec)')
    parser.add_argument('-dd', '--dr-dim', metavar='INT', required=False, default=20, type=int,
                        help='Number of dimensions to reduce to (default: 20)')
    parser.add_argument('-res', '--resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the the communities (default: 1.822)')
    parser.add_argument('-minr', '--min-resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the minimum size of the communities (default: 1.822)')
    parser.add_argument('-maxr', '--max-resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the maximum size of the communities (default: 1.822)')
    parser.add_argument('-ss', '--step-size', metavar='FLOAT', required=False, default=1.0, type=float,
                        help='Determines the step size to sweep resolution from min_resolution to max_resolution '
                             '(default: 1)')
    parser.add_argument('-nw', '--num-workers', metavar='INT', required=False, default=1, type=int,
                        help='Number of workers to run in parallel (default: 1, suggested: 25)')
    parser.add_argument('-nnm', '--num-neighbors-mi', metavar='INT', required=False, default=160, type=int,
                        help='Number of neighbors to build mutual information-based nearest neighbor graph '
                             '(default: 80)')
    parser.add_argument('-wl', '--walk-length', metavar='INT', required=False, default=15, type=int,
                        help='Length of random walks per source for graph embedding (default: 15)')
    parser.add_argument('-nl', '--num-walks', metavar='INT', required=False, default=20, type=int,
                        help='Number of random walks per source for graph embedding (default: 20)')
    parser.add_argument('-ws', '--window-size', metavar='INT', required=False, default=5, type=int,
                        help='Context window size of a random walk (default: 5)')
    parser.add_argument('-hp', '--hyper-p', metavar='FLOAT', required=False, default=1, type=float,
                        help='Hyperparameter p controls the likelihood of immediately traveling back to a node '
                             'recently traversed (default: 1)')
    parser.add_argument('-hq', '--hyper-q', metavar='FLOAT', required=False, default=1, type=float,
                        help='Hyperparameter q controls the likelihood of walking away from the previous '
                             'node (default: 1)')
    parser.add_argument('-nne', '--num-neighbors-eu', metavar='INT', required=False, default=25, type=int,
                        help='Number of neighbors to build euclidean distance-based nearest neighbor graph after '
                             'dimension reduction (default: 20)')
    parser.add_argument('-mc', '--meta-cell', required=False, action='store_true',
                        help='Create a MetaCell for each cell cluster.')
    
    parser.add_argument('-annef', '--ann-ef', metavar='INT', required=False, default=800, help='ef value of hnsw', type=int)
    parser.add_argument('-annm', '--ann-m', metavar='INT', required=False, default=4, help='M value of hnsw', type=int)
    parser.add_argument('-bpr', '--bin-power', metavar='INT', required=False, default=0, help='set the power index of the bin size for MI', type=int)
    parser.add_argument('-bsz', '--bin-size', metavar='INT', required=False, default=0, help='set the bin size for MI', type=int)

    parser.add_argument('-cldis', '--clustering-distance', metavar='STR', 
                        required=False, default='euclidean', help='euclidean/cosine')

    # parser.add_argument('-ha', '--harmony', metavar='FILE', required=False,
    #                     help='Path to a cell metadata file (tab-delimited text file) with "batch" as a column, '
    #                          'required for Harmony batch correction.')
    parser.add_argument('-sil', '--silhouette', metavar='INT', required=False, default=0, help='silhouette analysis(0) or not(other number)', type=int)
    parser.add_argument('-cs', '--consensus', metavar='STR', required=False, default='None', type=str,
                        choices=['None', 'CSPA', 'MCLA'], help='Consensus clustering methods. None means skip '
                                                               'consensus clustering;'
                                                               'CSPA is cluster-based similarity partitioning '
                                                               'algorithm; MCLA is meta-clustering algorithm. '
                                                               'Reference: '
                                                               'https://en.wikipedia.org/wiki/Consensus_clustering')
    return parser


def hnswlib_ann(data, num_neighbors=100, ef=200, M=16, num_jobs=1):
    dim = data.shape[1]
    num_elements = data.shape[0]
    ids = np.arange(num_elements)

    p = hnswlib.Index(space='mi', dim=dim)
    p.set_num_threads(num_jobs)
    p.init_index(max_elements = num_elements, ef_construction = ef, M = M)
    p.add_items(data, ids)
    p.set_ef(ef)
    labels, distances = p.knn_query(data, k = num_neighbors)
    return labels, distances



def mica_ge(args):
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    adata = pp.read_preprocessed_mat(args.input_file)
    frame = adata.to_df().astype(np.float32)
    del(adata)
    logging.info('(cells, genes): {}'.format(frame.shape))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building MI-based kNN graph on {} ...'.format(args.dr_modality))
    if args.dr_modality == 'gene':
        logging.info('Running HNSW ANN mode.')
        hnswlib.set_bin_power(args.bin_power)
        hnswlib.set_bin_size(args.bin_size)
        knn_indices, knn_dists = hnswlib_ann(frame.to_numpy(), 
                                             ef=args.ann_ef,
                                             M=args.ann_m,
                                             num_neighbors=args.num_neighbors_mi,
                                             num_jobs=args.num_workers)
            
    elif args.dr_modality == 'cell':
        logging.info('Running HNSW ANN mode.')
        hnswlib.set_bin_power(args.bin_power)
        hnswlib.set_bin_size(args.bin_size)
        knn_indices, knn_dists = hnswlib_ann(frame.T.to_numpy(), 
                                             ef=args.ann_ef,
                                             M=args.ann_m,
                                             num_neighbors=args.num_neighbors_mi, 
                                             num_jobs=args.num_workers)
            
    logging.info('kNN/ANN costs {}s'.format(time.time() - start))
    logging.info('kNN MI distance shape: {}'.format(knn_dists.shape))
    knn_graph = ng.general_graph_builder(knn_indices, knn_dists)
    logging.info('kNN graph number of nodes: {}'.format(len(knn_graph.nodes())))

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
    emb_file = '{}/knn_{}_graph_emb_on_{}_to_{}.txt'.format(args.output_dir, args.dr_method, args.dr_modality,
                                                            args.dr_dim)
    if args.dr_method == 'node2vec':
        if args.dr_modality == 'gene':
            dr.dim_reduce_node2vec_pecanpy(edgelist_file, emb_file, dim=args.dr_dim, num_jobs=args.num_workers,
                                           walk_len=args.walk_length, n_walks=args.num_walks,
                                           context_size=args.window_size, hyper_p=args.hyper_p, hyper_q=args.hyper_q)
            # wv = dr.dim_reduce_node2vec(knn_graph, dim=args.dr_dim, walk_len=10, n_walks=10)
            # print(wv)
        elif args.dr_modality == 'cell':
            dr.dim_reduce_node2vec_pecanpy(edgelist_file, emb_file, dim=args.dr_dim, num_jobs=args.num_workers,
                                           walk_len=args.walk_length, n_walks=args.num_walks,
                                           context_size=args.window_size, hyper_p=args.hyper_p, hyper_q=args.hyper_q)
    elif args.dr_method == 'deepwalk':
        # dr.dim_reduce_deepwalk(edgelist_file, emb_file, dim=args.dr_dim)
        sys.exit('Error - deepwalk has not been tested: {}'.format(args.dr_method))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))
    if args.dr_modality == 'cell':
        logging.info('All done.')
        sys.exit(0)

    start = time.time()
    logging.info('Performing clustering ...')
    mat_dr_df = pd.read_csv(emb_file, delimiter=' ', skiprows=1, index_col=0, names=np.arange(1, args.dr_dim+1))
    logging.info('(cells, genes): {}'.format(frame.shape))
    mat_dr_df.sort_index(inplace=True)
    mat_dr = mat_dr_df.to_numpy()
    logging.info('(cells, dimensions): {}'.format(mat_dr.shape))

    # Comment out this step as it does not work well based on a preliminary evaluation by Zhen Xie
    # if args.harmony:
    #     from harmony import harmonize
    #     start = time.time()
    #     logging.info('Performing harmony batch correction ...')
    #     cell_meta_df = pd.read_csv(args.harmony, sep='\t', index_col=0, header=0)
    #     mat_dr = harmonize(mat_dr, cell_meta_df, batch_key='batch')
    #     end = time.time()
    #     runtime = end - start
    #     logging.info('Done. Runtime: {} seconds'.format(runtime))
    logging.info('Building Neighbor Graph...')
    G = ng.build_graph(mat_dr, 
                       dis_metric=args.clustering_distance, 
                       num_neighbors=args.num_neighbors_eu, 
                       num_jobs=args.num_workers)
    
    edgelist_file = '{}/sklearn_knn_graph_edgelist.txt'.format(args.output_dir)
    nx.write_edgelist(G, edgelist_file)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))
    start = time.time()
    logging.info('Louvain clustering...')
    if (args.max_resolution == args.min_resolution) and (args.resolution != args.min_resolution):
        partition_resolutions = cl.graph_clustering_parallel(G, min_resolution=args.resolution,
                                                             max_resolution=args.resolution,
                                                             step_size=args.step_size, num_workers=args.num_workers)
    else:
        partition_resolutions = cl.graph_clustering_parallel(G, min_resolution=args.min_resolution,
                                                             max_resolution=args.max_resolution,
                                                             step_size=args.step_size, num_workers=args.num_workers)
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    # Create MetaCells
    if args.meta_cell:
        start = time.time()
        logging.info('Creating Metacells ...')
        mc.create_metacells(partition_resolutions, G, frame.iloc[mat_dr_df.index,:].index, args.output_dir)
        end = time.time()
        runtime = end - start
        logging.info('Done. Runtime: {} seconds'.format(runtime))

    aggs = []
    if args.consensus != 'None':
        start = time.time()
        logging.info('Group partitions based on the number of clusters and consensus clustering ...')
        partitions = [par for (par, reso) in partition_resolutions]
        cluster_dict = cs.group_partition(partitions, frame.index)
        for num_cluster in cluster_dict.keys():
            agg, out_f = cs.consensus_sc3(cluster_dict[num_cluster], num_cluster, 'consensus')
            aggs.append((agg, num_cluster))
        end = time.time()
        runtime = end - start
        logging.info('Done. Runtime: {} seconds'.format(runtime))
    else:
        for par, reso in partition_resolutions:
            labels = [x + 1 for x in par.values()]
            clustering_res = pd.DataFrame(data=labels, index=frame.iloc[mat_dr_df.index,:].index, columns=["label"])
            aggs.append((clustering_res, reso))

    logging.info('Visualizing clustering results using {}'.format(args.visual_method))
    aggs_embed = []
    for partition, num_cluster_metric in aggs:      # num_cluster_metric is either resolution or num_cluster
        if args.consensus == 'None':
            num_cluster_metric = round(num_cluster_metric, 5)
        agg_embed = vs.visual_embed(partition, mat_dr, args.output_dir, suffix=num_cluster_metric,
                                    visual_method=args.visual_method, num_works=args.num_workers,
                                    min_dist=args.min_dist)
        aggs_embed.append((agg_embed, num_cluster_metric))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))


    if args.silhouette != 0:
        start = time.time()
        logging.info('Optimal number of clusters analysis ...')
        with open('{}/silhouette_avg.txt'.format(args.output_dir), 'w') as fout:
            if args.consensus != 'None':
                fout.write('dimension\tnum_clusters\tsilhouette_avg\n')
                for agg, num_clusters in aggs_embed:
                    logging.info('number of clusters: {}'.format(num_clusters))
                    logging.info(agg)
                    silhouette_avg = cl.silhouette_analysis(agg, num_clusters, agg.loc[:, ['X', 'Y']], args.output_dir)
                    fout.write('{}\t{}\t{}\n'.format(args.dr_dim, num_clusters, silhouette_avg))
            else:
                fout.write('dimension\tresolution\tnum_clusters\tsilhouette_avg\n')
                for agg, resolution in aggs_embed:
                    resolution = round(resolution, 5)
                    logging.info('resolution: {}'.format(resolution))
                    # logging.info(agg)
                    num_clusters = len(set(agg['label']))
                    silhouette_avg = cl.silhouette_analysis(agg, num_clusters, agg.loc[:, ['X', 'Y']], args.output_dir,
                                                            resolution=resolution)
                    fout.write('{}\t{}\t{}\t{}\n'.format(args.dr_dim, resolution, num_clusters, silhouette_avg))
        end = time.time()
        runtime = end - start
        logging.info('Done. Runtime: {} seconds'.format(runtime))

    # start = time.time()
    # out_h5_file = '{}/clustered.h5ad'.format(args.output_dir)
    # logging.info('Write h5ad file {} ...'.format(out_h5_file))
    # pp.write_h5(adata, out_h5_file)
    # end = time.time()
    # runtime = end - start
    # logging.info('Done. Runtime: {} seconds'.format(runtime))

    logging.info('All done.')


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# MICA GE version uses a graph embedding method for dimension reduction on MI-kNN graph.

import sys
import logging
import time
import argparse
import pickle
import numpy as np
import pandas as pd
import pathlib
from sklearn.cluster import KMeans
from MICA.lib import transform as trans
from MICA.lib import mutual_info as mi
from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs
from MICA.lib import consensus as cs

def main():
    head_description = '''MICA - Mutual Information-based Clustering Analysis tool. This version uses a 
    multidimensional scaling method for dimension reduction.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                                     description=head_description)
    
    parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                        help='Path to an input file (h5ad file or tab-delimited text file)')
    parser.add_argument('-o', '--output-dir', metavar='DIR', required=True, help='Path to final output directory')
    parser.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                        help='Visualization method UMAP or t-SNE (default: UMAP)')
    parser.add_argument('-md', '--min-dist', metavar='FLOAT', required=False, default=0.6, type=float,
                        help='min_dist parameter in UMAP, minimum distance of points in the embedded space '
                             '(default: 0.6)')

    parser = add_mds_arguments(parser)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mica_mds(args)


def add_mds_arguments(parser):
    parser.add_argument('-dm', '--dr-method', metavar='STR', default='MDS', required=False,
                        help='Transformation method used for dimension reduction '
                             '[MDS | PCA] (default: MDS)')
    parser.add_argument('-dd', '--dr-dim', metavar='INT', required=False, default=19, type=int,
                        help='Number of dimensions to reduce to (default: 20)')
    parser.add_argument('-nw', '--num-workers', metavar='INT', required=False, default=1, type=int,
                        help='Number of workers to run in parallel (default: 1, suggested: 25)')

    # parser.add_argument('-nc', '--num-clusters', metavar='INT', nargs='+', required=False, default=0, type=int,
    #                       help='Number of clusters to be specified in kmeans')
    parser.add_argument('-nck', '--num-clusters-k', metavar='INT', default=4, required=False, type=int,
                                 help='Number of clusters to be specified in kmeans')
    parser.add_argument('-le', '--louvain-enable', metavar='INT', required=False, default=0, help='enable knn-louvain clustering or not(0)', type=int)
    parser.add_argument('-nn', '--num-neighbors', metavar='INT', required=False, default=20, type=int,
                        help='Number of neighbors to build euclidean distance-based nearest neighbor graph after '
                            'dimension reduction (default: 20)')
    parser.add_argument('-cldis', '--clustering-distance', metavar='STR', 
                        required=False, default='euclidean', help='euclidean/cosine')
    parser.add_argument('-res', '--resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the the communities (default: 1.822)')
    parser.add_argument('-minr', '--min-resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the minimum size of the communities (default: 1.822)')
    parser.add_argument('-maxr', '--max-resolution', metavar='FLOAT', required=False, default=1.822, type=float,
                        help='Determines the maximum size of the communities (default: 1.822)')
    parser.add_argument('-ss', '--step-size', metavar='FLOAT', required=False, default=1, type=float,
                        help='Determines the step size to sweep resolution from min_resolution to max_resolution '
                             '(default: 1)')

    parser.add_argument('-bpr', '--bin-power', metavar='INT', required=False, default=0, help='set the power index of the bin size for MI', type=int)
    parser.add_argument('-bsz', '--bin-size', metavar='INT', required=False, default=0, help='set the bin size for MI', type=int)


    parser.add_argument('-sil', '--silhouette', metavar='INT', required=False, default=0, help='silhouette analysis(other numbers) or not(0)', type=int)
    parser.add_argument('-cs', '--consensus', metavar='STR', required=False, default='None', type=str,
                        choices=['None', 'CSPA', 'MCLA'], help='Consensus clustering methods. None means skip '
                                                               'consensus clustering;'
                                                               'CSPA is cluster-based similarity partitioning '
                                                               'algorithm; MCLA is meta-clustering algorithm. '
                                                               'Reference: '
                                                               'https://en.wikipedia.org/wiki/Consensus_clustering')
    return parser


def mica_mds(args):
    start = time.time()
    logging.info('Read preprocessed expression matrix ...')
    adata = pp.read_preprocessed_mat(args.input_file)
    frame = adata.to_df()
    data = frame.to_numpy()
    logging.info('(cells, genes): {}'.format(frame.shape))
    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))

    start = time.time()
    logging.info('Building MI matrix ...')
    if args.bin_size != 0:
        bins = args.bin_size
    elif args.bin_power != 0:
        bins = int(data.shape[1] ** (1/args.bin_power))
    else:
        bins = int(data.shape[1] ** (1/3.0))
    mi_matrix = mi.mi_norm(mi.mi_calc(data, bins))
    
    logging.info('MI calculation costs {}s'.format(time.time() - start))
    if args.dr_method == 'MDS':
        logging.info('Performing MDS transformation ...')
        mi_mds = trans.mds_trans(mi_matrix)[:, :args.dr_dim]
    elif args.dr_method == 'PCA':
        logging.info('Performing PCA transformation ...')
        mi_mds = trans.pca_trans(mi_matrix, args.dr_dim)
    else:
        logging.info('Performing MDS transformation ...')
        mi_mds = trans.mds_trans(mi_matrix)[:, :args.dr_dim]
    logging.info('dim-reduced data shape {}'.format(mi_mds.shape))

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    with open(args.output_dir + '/mi_reduced.pkl', 'wb') as f:
        pickle.dump(mi_mds, f)
        f.close()
    with open(args.output_dir + '/mi_info.pkl', 'wb') as f:
        pickle.dump(mi_matrix, f)
        f.close()

    end = time.time()
    runtime = end - start
    logging.info('Done. Runtime: {} seconds'.format(runtime))


    start = time.time()
    logging.info('Performing clustering ...')
    logging.info('(cells, genes): {}'.format(frame.shape))

    if not args.louvain_enable:
        logging.info('Performing Kmeans clustering for # of {}...'.format(args.num_clusters_k))
        kms = KMeans(n_clusters=args.num_clusters_k)
        kms.fit(mi_mds)
        partition = pd.DataFrame(data=[i + 1 for i in list(kms.labels_)], index=frame.index, columns=["label"])

        agg_embed = vs.visual_embed(partition, mi_mds,
                                    args.output_dir, 
                                    suffix=args.num_clusters_k,
                                    visual_method=args.visual_method, 
                                    num_works=args.num_workers,
                                    min_dist=args.min_dist)
        end = time.time()
        runtime = end - start
        logging.info('Done. Runtime: {} seconds'.format(runtime))
    
    else:
        logging.info('Performing KNN clustering ...')
        G = ng.build_graph(mi_mds, 
                           dis_metric=args.clustering_distance, 
                           num_neighbors=args.num_neighbors, 
                           num_jobs=args.num_workers)
        stamp1 = time.time()
        logging.info('Graph building cost {} seconds'.format(stamp1 - start))
        if (args.max_resolution == args.min_resolution) and (args.resolution != args.min_resolution):
            partition_resolutions = cl.graph_clustering_parallel(G, min_resolution=args.resolution,
                                                                    max_resolution=args.resolution,
                                                                    step_size=args.step_size, 
                                                                    num_workers=args.num_workers)
        else:
            partition_resolutions = cl.graph_clustering_parallel(G, min_resolution=args.min_resolution,
                                                                    max_resolution=args.max_resolution,
                                                                    step_size=args.step_size, 
                                                                    num_workers=args.num_workers)
        end = time.time()
        logging.info('Louvain clustering cost {} seconds'.format(end - stamp1))
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
                clustering_res = pd.DataFrame(data=labels, index=frame.index, columns=["label"])
                aggs.append((clustering_res, reso))

        logging.info('Visualizing clustering results using {}'.format(args.visual_method))
        aggs_embed = []
        for partition, num_cluster_metric in aggs:      # num_cluster_metric is either resolution or num_cluster
            if args.consensus == 'None':
                num_cluster_metric = round(num_cluster_metric, 5)
            agg_embed = vs.visual_embed(partition, mi_mds, args.output_dir, suffix=num_cluster_metric,
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

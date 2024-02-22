#!/usr/bin/env python3
import glob
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score
import seaborn as sns


def create_summary_table():
    directory = '/Users/lding/Documents/MICA/tests/SilverStd/PBMC_20k/nne'
    pathname = directory + "/nne_MICA_GE_out/**/MICA_GE.out"
    files = glob.glob(pathname, recursive=True)
    true_label_file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_true_label.txt'
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    summary = []
    for f in files:
        tokens = f.split('/')
        test_run_dirname = tokens[10]
        print(test_run_dirname)
        paras_tokens = test_run_dirname.split('_')
        nn = int(paras_tokens[1])
        nj = int(paras_tokens[3])
        # Parse stdout file
        with open(f, 'r') as mica_out:
            for line in mica_out:
                if line.find('Max Memory') != -1:
                    line_tokens = line.split()
                    max_mem = int(line_tokens[3])
                elif line.find('Run time') != -1:
                    line_tokens = line.split()
                    run_time = int(line_tokens[3])
                    break

        # Parse silhouette_avg file
        silhouette_avg_file = directory + "/nne_silhouette_avg/" + test_run_dirname + "/silhouette_avg.txt"
        resolution = None
        with open(silhouette_avg_file, 'r') as sil_f:
            for line in sil_f:
                line_tokens = line.split()
                if line_tokens[2] in ['9', '10', '11']:
                    resolution = float(line_tokens[1])
                    num_clusters = int(line_tokens[2])
                    silhouette = float(line_tokens[3])

                    # Parse clustering_UMAP file
                    cluster_umap_file = '{}/nne_clustering_UMAP/{}/clustering_UMAP_euclidean_20_{}.txt'.format(
                                        directory, test_run_dirname, resolution)
                    predict_label = pd.read_csv(cluster_umap_file, delimiter='\t', index_col=0)
                    merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                    ari = adjusted_rand_score(merged['label_x'], merged['label_y'])

                    summary.append((nn, nj, max_mem, run_time, resolution, num_clusters, silhouette, ari))
        if resolution is None:
            continue

    summary_df = pd.DataFrame(data=summary, columns=['num_neighbor_euclidean', 'num_workers', 'max_mem', 'run_time',
                                                     'resolution', 'num_clusters', 'silhouette', 'ARI'])
    # print(summary_df)
    summary_df.to_csv('{}/summary.txt'.format(directory), sep='\t')
    return


def create_summary_table_hpc():
    directory = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/' \
                'computing_performance/pecanpy/hyperparameter'
    pathname = directory + "/**/MICA_GE.out"
    files = glob.glob(pathname, recursive=True)
    true_label_file = '/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/' \
                      'PBMC_20k/PBMC_20k_true_label.txt'
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    summary = []
    for f in files:
        tokens = f.split('/')
        test_run_dirname = tokens[14]
        print(test_run_dirname)
        paras_tokens = test_run_dirname.split('_')
        phyper = float(paras_tokens[1])
        qhyper = float(paras_tokens[3])
        nw = int(paras_tokens[5])
        # Parse stdout file
        with open(f, 'r') as mica_out:
            for line in mica_out:
                if line.find('Max Memory') != -1:
                    line_tokens = line.split()
                    max_mem = int(line_tokens[3])
                elif line.find('Run time') != -1:
                    line_tokens = line.split()
                    run_time = int(line_tokens[3])
                    break

        # Parse silhouette_avg file
        silhouette_avg_file = directory + '/' + test_run_dirname + "/silhouette_avg.txt"
        resolution = None
        with open(silhouette_avg_file, 'r') as sil_f:
            for line in sil_f:
                line_tokens = line.split()
                if line_tokens[2] in ['9', '10', '11']:
                    resolution = float(line_tokens[1])
                    num_clusters = int(line_tokens[2])
                    silhouette = float(line_tokens[3])

                    # Parse clustering_UMAP file
                    cluster_umap_file = '{}/{}/clustering_UMAP_euclidean_20_{}.txt'.format(
                                        directory, test_run_dirname, resolution)
                    predict_label = pd.read_csv(cluster_umap_file, delimiter='\t', index_col=0)
                    merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                    ari = adjusted_rand_score(merged['label_x'], merged['label_y'])

                    # print((wlen, nwalk, ws, max_mem, run_time, resolution, num_clusters, silhouette, ari))
                    summary.append((phyper, qhyper, nw, max_mem, run_time, resolution, num_clusters, silhouette, ari))
        if resolution is None:
            continue

    summary_df = pd.DataFrame(data=summary, columns=['hyperparameter_p', 'hyperparameter_q', 'num_workers', 'max_mem',
                                                     'run_time', 'resolution', 'num_clusters', 'silhouette', 'ARI'])
    # print(summary_df)
    summary_df.to_csv('{}/summary.txt'.format(directory), sep='\t')
    return


if __name__ == "__main__":
    # create_summary_table()
    create_summary_table_hpc()

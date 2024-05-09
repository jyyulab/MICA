#!/usr/bin/env python3

import os
import pathlib
import pandas as pd
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


def create_cmds_ge(root_dir, level, author, input_file, max_resolution):
    test_dir = '{}/tests/{}/{}/ge'.format(root_dir, level, author)
    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)
    output_dir = '{}/outputs/{}/{}/ge'.format(root_dir, level, author)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    cmd_file = '{}/cmd_eval.sh'.format(test_dir)
    input_file = '{}/datasets/{}/{}/{}'.format(root_dir, level, author, input_file)
    with open(cmd_file, 'w') as fout:
        for dim in range(8, 100, 4):
            out_dir = '{}/dim_{}'.format(output_dir, dim)
            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)
            cmd = 'mica ge -i {} -o {} -dd {} -ar {}'.format(input_file, out_dir, dim, max_resolution)
            fout.write(cmd)
            fout.write('\n')
    print('cmd file: {}'.format(cmd_file))
    print('Done')


def create_cmds_mds(root_dir, level, author, input_file, num_clusters):
    test_dir = '{}/tests/{}/{}/mds'.format(root_dir, level, author)
    if not os.path.isdir(test_dir):
        pathlib.Path(test_dir).mkdir(parents=True, exist_ok=True)
    output_dir = '{}/outputs/{}/{}/mds'.format(root_dir, level, author)
    if not os.path.isdir(output_dir):
        pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    cmd_file = '{}/cmd_eval.sh'.format(test_dir)
    input_file = '{}/datasets/{}/{}/{}'.format(root_dir, level, author, input_file)

    with open(cmd_file, 'w') as fout:
        cmd = 'mica mds -i {} -o {} -pn {} -nc {} -dk 12 13 14 15 16 17 18 19\n'.format(input_file, output_dir, author,
                                                                                        num_clusters)
        print(cmd)
        fout.write(cmd)
    print('cmd file: {}'.format(cmd_file))
    print('Done')


def create_cmds_auto(root_dir, level, author, input_file):
    pass


def calc_ARIs_ge(root_dir, level, author, num_clusters):
    output_dir = '{}/outputs/{}/{}/ge_1_3'.format(root_dir, level, author)
    # output_dir = '{}/outputs/{}/{}/old'.format(root_dir, level, author)
    summary_file = '{}/summary_ge.txt'.format(output_dir)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
    print(true_label_file)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    best_ari = 0.0
    with open(summary_file, 'w') as fout:
        for dim in range(8, 100, 4):
            clustering_out_dir = '{}/dim_{}'.format(output_dir, dim)
            # for reso in np.arange(0.4, 10.1, 0.4):
            for reso in np.arange(0.2, 10.0, 0.4):
                reso_round = np.round(reso, 2)
                # clustering_out_dir = '{}/{}_{}'.format(output_dir, dim, reso_round)
                predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(clustering_out_dir, dim,
                                                                                     reso_round)
                # predict_label_file = '{}/clustering_umap_euclidean_{}_{}.txt'.format(clustering_out_dir, dim,
                #                                                                      reso_round)
                predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
                predict_num_clusters = len(set(predict_label['label']))
                # print(predict_num_clusters)
                # print('dim_{}_reso_{}_numCluster_{}'.format(dim, reso_round, predict_num_clusters))
                if predict_num_clusters != num_clusters:
                   continue
                # new_index = [int(s.replace('V', '')) for s in predict_label.index]
                # predict_label.index = new_index
                merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                # print(merged)
                # ari = adjusted_rand_score(merged['cell'], merged['label_y'])
                ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
                print('{} {} {} {}'.format(predict_num_clusters, dim, reso, ari))
                if ari > best_ari:
                    best_ari = ari
                    best_ari_str = '{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ari)
                fout.write('{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ari))
                # break
            # break
    print('GE best ARI: {}'.format(best_ari_str))
    print('Done')


def calc_ARIs_mds(root_dir, level, author, num_clusters):
    output_dir = '{}/outputs/{}/{}'.format(root_dir, level, author)
    summary_file = '{}/summary_mds.txt'.format(output_dir)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    print(true_label)
    with open(summary_file, 'w') as fout:
        predict_label_file = '{}/mds_1_3/{}_scDHA_k{}_umap_ClusterMem.txt'.format(output_dir, author, num_clusters)
        print(predict_label_file)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
        # print(predict_label)
        # new_index = [int(s.replace('V', '')) for s in predict_label.index]
        # predict_label.index = new_index
        # merged = true_label.merge(predict_label, left_on='cell', right_index=True)
        merged = true_label.merge(predict_label, left_on='cell', right_index=True)
        print(merged)
        # ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
        ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
        head_line = 'author\tnum_clusters\tARI\n'
        print('MDS')
        print(head_line, end='')
        fout.write(head_line)
        out_line = '{}\t{}\t{}\n'.format(author, num_clusters, ari)
        print(out_line, end='')
        fout.write(out_line)
    print('Done')


def calc_AMIs_mds(root_dir, level, author, num_clusters):
    output_dir = '{}/outputs/{}/{}'.format(root_dir, level, author)
    summary_file = '{}/summary_mds.txt'.format(output_dir)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
    # true_label_file = '{}/datasets/{}/{}/scGNN/{}_cell_label.csv'.format(root_dir, level, author, author)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    print(true_label)
    with open(summary_file, 'w') as fout:
        predict_label_file = '{}/mds_1_3_silhouette/{}_k{}_umap_ClusterMem.txt'.format(output_dir, author, num_clusters)
        print(predict_label_file)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
        # print(predict_label)
        # new_index = [int(s.replace('V', '')) for s in predict_label.index]
        # predict_label.index = new_index
        # merged = true_label.merge(predict_label, left_on='cell', right_index=True)
        merged = true_label.merge(predict_label, left_index=True, right_index=True)
        print(merged)
        # print(set(merged['label_x']))
        # ari = adjusted_rand_score(merged['label_x'], merged['label_y'])
        ami = adjusted_mutual_info_score(merged['label_x'], merged['label_y'])
        head_line = 'author\tnum_clusters\tAMI\n'
        print('MDS')
        print(head_line, end='')
        fout.write(head_line)
        out_line = '{}\t{}\t{}\n'.format(author, num_clusters, ami)
        print(out_line, end='')
        fout.write(out_line)
    print('Done')


def calc_AMIs_ge(root_dir, level, author, num_clusters):
    output_dir = '{}/outputs/{}/{}/ge'.format(root_dir, level, author)
    # output_dir = '{}/outputs/{}/{}/old'.format(root_dir, level, author)
    summary_file = '{}/summary_ge.txt'.format(output_dir)
    if os.path.isfile(summary_file):
        os.remove(summary_file)

    true_label_file = '{}/datasets/{}/{}/{}_true_label.txt'.format(root_dir, level, author, author)
    print(true_label_file)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    best_ami = 0.0
    with open(summary_file, 'w') as fout:
        for dim in range(8, 100, 4):
            clustering_out_dir = '{}/dim_{}'.format(output_dir, dim)
            # for reso in np.arange(0.4, 10.1, 0.4):
            for reso in np.arange(0.2, 10.0, 0.4):
                reso_round = np.round(reso, 2)
                # clustering_out_dir = '{}/{}_{}'.format(output_dir, dim, reso_round)
                predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(clustering_out_dir, dim,
                                                                                     reso_round)
                # predict_label_file = '{}/clustering_umap_euclidean_{}_{}.txt'.format(clustering_out_dir, dim,
                #                                                                      reso_round)
                predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
                predict_num_clusters = len(set(predict_label['label']))
                # print(predict_num_clusters)
                # print('dim_{}_reso_{}_numCluster_{}'.format(dim, reso_round, predict_num_clusters))
                if predict_num_clusters != num_clusters:
                   continue
                # new_index = [int(s.replace('V', '')) for s in predict_label.index]
                # predict_label.index = new_index
                merged = true_label.merge(predict_label, left_on='cell', right_index=True)
                # print(merged)
                # print(set(merged['class_label']))
                # print(len(set(merged['class_label'])))
                # print(set(merged['label']))
                # print(len(set(merged['label'])))
                # ami = adjusted_mutual_info_score(merged['cell'], merged['label_y'])
                # ami = adjusted_mutual_info_score(merged['label_x'], merged['label_y'])
                ami = adjusted_mutual_info_score(merged['subclass_label'], merged['label'])
                print('{} {} {} {}'.format(predict_num_clusters, dim, reso, ami))
                if ami > best_ami:
                    best_ami = ami
                    best_ami_str = '{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ami)
                fout.write('{}\t{}\t{}\t{}\n'.format(author, dim, reso_round, ami))
                # break
            # break
    print('GE best ARI: {}'.format(best_ami_str))
    print('Done')

#!/usr/bin/env python3

from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
from MICA.lib.aux_utils import run_shell_command
import numpy as np
import os


def GSE71585_Tasic(dim, reso):
    true_label_file = '{}/datasets/GSE71585_Tasic/GSE71585_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['type'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['type']]
    true_label['label'] = labels
    gse71582 = '{}/datasets/GSE71585_Tasic/GSE71585_preprocessed.h5ad'.format(mica_data_path)
    gse71582_out = '{}/outputs/GSE71585_Tasic/new'.format(mica_data_path)
    if not os.path.isdir(gse71582_out):
        os.makedirs(gse71582_out)
    # cmd = 'mica -i {} -o {} -d {} -e {} -s 0.2'.format(gse71582, gse71582_out, dim, reso)
    # print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(gse71582_out, dim, reso)
    predict_label = pd.read_csv(predict_label_file, delimiter='\t')
    return adjusted_rand_score(true_label['label'], predict_label['label'])


def GSE60361_Ziesel(dim, reso):
    true_label_file = '{}/datasets/GSE60361_Ziesel/Ziesel_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['type'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['type']]
    true_label['label'] = labels
    ziesel = '{}/datasets/GSE60361_Ziesel/Ziesel_MICA_input.txt'.format(mica_data_path)
    Ziesel_out = '{}/outputs/GSE60361_Ziesel/new'.format(mica_data_path)
    if not os.path.isdir(Ziesel_out):
        os.makedirs(Ziesel_out)
    # cmd = 'mica -i {} -o {} -d {} -e {} -s 0.2'.format(ziesel, Ziesel_out, dim, reso)
    # print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(Ziesel_out, dim, reso)
    predict_label = pd.read_csv(predict_label_file, delimiter='\t')
    return adjusted_rand_score(true_label['label'], predict_label['label'])


def GSE75688_Chung(dim, reso):
    true_label_file = '{}/datasets/GSE75688_Chung/GSE75688_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['type'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['type']]
    true_label['label'] = labels
    Chung = '{}/datasets/GSE75688_Chung/GSE75688_Chung_preprocessed.h5ad'.format(mica_data_path)
    Chung_out = '{}/outputs/GSE75688_Chung/new'.format(mica_data_path)
    if not os.path.isdir(Chung_out):
        os.makedirs(Chung_out)
    # cmd = 'mica -i {} -o {} -d {} -e {} -s 0.2'.format(Chung, Chung_out, dim, reso)
    # print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(Chung_out, dim, reso)
    predict_label = pd.read_csv(predict_label_file, delimiter='\t')
    return adjusted_rand_score(true_label['label'], predict_label['label'])


def PBMC20k_MDS():
    for dim in range(9, 14):
        true_label_file = '{}/datasets/PBMC_20k/PBMC_sorted_20K_true_label.txt'.format(mica_data_path)
        true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
        cell_type_label_dict = dict()
        for i, v in enumerate(set(true_label['type'])):
            cell_type_label_dict[v] = i
        labels = [cell_type_label_dict[ct] for ct in true_label['type']]
        true_label['label'] = labels
        # PBMC20k_out = '{}/outputs/PBMC_20k_MDS/mica-9a7645e0-5db0-4e37-8a2a-1a68d3471395'.format(mica_data_path)
        PBMC20k_out = '{}/outputs/PBMC_20k_MDS/mica-ba57c9e0-a2a6-44a0-a445-d46ec51a21d8'.format(mica_data_path)
        predict_label_file = '{}/cwl_lsf_k{}_tsne_ClusterMem.txt'.format(PBMC20k_out, dim)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t')
        print(dim)
        print(adjusted_rand_score(true_label['label'], predict_label['label']))


def GSE71585_Tasic_MDS():
    for dim in range(2, 11):
        true_label_file = '{}/datasets/GSE71585_Tasic/GSE71585_true_label.txt'.format(mica_data_path)
        true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
        cell_type_label_dict = dict()
        for i, v in enumerate(set(true_label['type'])):
            cell_type_label_dict[v] = i
        labels = [cell_type_label_dict[ct] for ct in true_label['type']]
        true_label['label'] = labels
        Tasic_out = '{}/outputs/GSE71585_Tasic/MDS/mica-bada1315-e689-4858-b2b3-c773bd3e95c5'.format(mica_data_path)
        predict_label_file = '{}/cwl_lsf_k{}_tsne_ClusterMem.txt'.format(Tasic_out, dim)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t')
        print(dim)
        # print('{}'.format(adjusted_rand_score(true_label['label'], predict_label['label'])))


def GSE60361_Ziesel_MDS():
    for dim in range(5, 11):
        true_label_file = '{}/datasets/GSE60361_Ziesel/Ziesel_true_label.txt'.format(mica_data_path)
        true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
        cell_type_label_dict = dict()
        for i, v in enumerate(set(true_label['type'])):
            cell_type_label_dict[v] = i
        labels = [cell_type_label_dict[ct] for ct in true_label['type']]
        true_label['label'] = labels
        Ziesel_out = '{}/outputs/GSE60361_Ziesel/MDS/mica-45a2fb5c-ce90-4472-aada-4de5097623c1'.format(mica_data_path)
        predict_label_file = '{}/cwl_lsf_k{}_tsne_ClusterMem.txt'.format(Ziesel_out, dim)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t')
        print('{}\t{}'.format(dim, adjusted_rand_score(true_label['label'], predict_label['label'])))


def GSE75688_Chung_MDS():
    for dim in range(3, 8):
        true_label_file = '{}/datasets/GSE75688_Chung/GSE75688_true_label.txt'.format(mica_data_path)
        true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
        cell_type_label_dict = dict()
        for i, v in enumerate(set(true_label['type'])):
            cell_type_label_dict[v] = i
        labels = [cell_type_label_dict[ct] for ct in true_label['type']]
        true_label['label'] = labels
        Chung_out = '{}/outputs/GSE75688_Chung/MDS/mica-48447539-182e-4d63-831e-08a664e5cded'.format(mica_data_path)
        predict_label_file = '{}/cwl_lsf_k{}_tsne_ClusterMem.txt'.format(Chung_out, dim)
        predict_label = pd.read_csv(predict_label_file, delimiter='\t')
        print('{}\t{}'.format(dim, adjusted_rand_score(true_label['label'], predict_label['label'])))


def Buttner(dim, reso):
    true_label_file = '{}/datasets/GoldernStd/buettner/buettner_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['label'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['label']]
    true_label['label'] = labels
    # print(true_label)
    yan = '{}/datasets/GoldernStd/buettner/Buttner_MICA_input.txt'.format(mica_data_path)
    yan_out = '{}/outputs/GoldernStd/Buttner/MICA_GE'.format(mica_data_path)

    if not os.path.isdir(yan_out):
        os.makedirs(yan_out)
    cmd = 'mica -i {} -o {} -d {} -e {} -s 0.2'.format(yan, yan_out, dim, reso)
    print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    # predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(yan_out, dim, reso)
    # predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
    # print(predict_label)
    # return adjusted_rand_score(true_label['label'], predict_label['label'])


def Kolod(dim, reso):
    true_label_file = '{}/datasets/GoldernStd/kolod/kolod_true_label.txt'.format(mica_data_path)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    cell_type_label_dict = dict()
    for i, v in enumerate(set(true_label['label'])):
        cell_type_label_dict[v] = i
    labels = [cell_type_label_dict[ct] for ct in true_label['label']]
    true_label['label'] = labels
    # print(true_label)
    yan = '{}/datasets/GoldernStd/kolod/Kolod_MICA_input.txt'.format(mica_data_path)
    yan_out = '{}/outputs/GoldernStd/Kolod/MICA_GE'.format(mica_data_path)

    if not os.path.isdir(yan_out):
        os.makedirs(yan_out)
    cmd = 'mica -i {} -o {} -d {} -e {} -s 0.2'.format(yan, yan_out, dim, reso)
    print(cmd)
    # run_shell_command(cmd)
    # print('Done')
    # predict_label_file = '{}/clustering_UMAP_euclidean_{}_{}.txt'.format(yan_out, dim, reso)
    # predict_label = pd.read_csv(predict_label_file, delimiter='\t', index_col=0)
    # print(predict_label)
    # return adjusted_rand_score(true_label['label'], predict_label['label'])


def create_cmds():
    for dim in range(8, 100, 4):
        for reso in np.arange(0.2, 10.0, 0.4):
            reso_round = np.round(reso, 1)
            # GSE71585_Tasic(dim, 10.0)
            # GSE60361_Ziesel(dim, 10.0)
            # GSE75688_Chung(dim, 10.0)
            # PBMC20k(reso_round)
            # Yan(dim, reso_round)
            # Buttner(dim, reso_round)


def create_cmds_NNDescent():
    for dim in range(8, 100, 4):
        for reso in np.arange(0.2, 10.0, 0.4):
            reso_round = np.round(reso, 1)
            # GSE71585_Tasic(dim, reso_round)
            GSE60361_Ziesel(dim, reso_round)
            # GSE75688_Chung(dim, reso_round)
            # PBMC20k(dim, reso_round)
            # PBMC20k(reso_round)


def summary():
    summary_file = '{}/outputs/summary.txt'.format(mica_data_path)
    if os.path.isfile(summary_file):
        os.remove(summary_file)
    with open('{}/outputs/summary.txt'.format(mica_data_path), 'w') as fout:
        for dim in range(8, 100, 4):
            for reso in np.arange(1.0, 10.0, 0.4):
                reso_round = np.round(reso, 1)
                print('DR dimension: {}\t louvain: {}\n'.format(dim, reso_round))
                # ari = GSE71585_Tasic(dim, reso_round)
                # print('ari: {}\n'.format(ari))
                # fout.write('GSE71585\t{}\t{}\t{}\n'.format(dim, reso_round, ari))
                #
                # ari = GSE60361_Ziesel(dim, reso_round)
                # print('ari: {}\n'.format(ari))
                # fout.write('Ziesel\t{}\t{}\t{}\n'.format(dim, reso_round, ari))
                #
                # ari = GSE75688_Chung(dim, reso_round)
                # print('ari: {}\n'.format(ari))
                # fout.write('Chung\t{}\t{}\t{}\n'.format(dim, reso_round, ari))

                # ari = Yan(dim, reso_round)
                # print('ari: {}\n'.format(ari))
                # fout.write('Yan\t{}\t{}\t{}\n'.format(dim, reso_round, ari))
                break
            break


def summary_no_ge():
    summary_file = '{}/outputs/summary_no_ge.txt'.format(mica_data_path)
    if os.path.isfile(summary_file):
        os.remove(summary_file)
    with open('{}/outputs/summary_no_ge.txt'.format(mica_data_path), 'w') as fout:
        for reso in np.arange(0.2, 10.0, 0.4):
            reso_round = np.round(reso, 1)

            # break


if __name__ == "__main__":
    mica_data_path = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    create_cmds()
    # summary()
    # summary_no_ge()
    # PBMC20k_MDS()
    # GSE71585_Tasic_MDS()
    # GSE60361_Ziesel_MDS()
    # GSE75688_Chung_MDS()

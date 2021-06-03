#!/usr/bin/env python3

import os
import glob
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score


def run_MDS(author):
    cluster_num = 5
    root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    test_dir = '{}/tests/SilverStd/{}/MICA_MDS/cmd_files'.format(root_dir, author)
    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)
    # input_file = '{}/datasets/SilverStd/{}/{}_MICA_input.txt'.format(root_dir, author, author)
    input_file = '{}/datasets/SilverStd/{}/{}_preprocessed.h5ad'.format(root_dir, author, author)
    for dim in range(1, 51):
        output_dir = '{}/outputs/SilverStd/{}/MICA_MDS/dim_{}'.format(root_dir, author, dim)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        with open('{}/test_MICA_MDS_lsf_{}.sh'.format(test_dir, dim), 'w') as fout:
            fout.write('#BSUB -q compbio\n#BSUB -P {} \n#BSUB -oo {}_{}.out -eo {}_{}.err\n'.format(author, author,
                                                                                                    dim, author, dim))
            fout.write('mica lsf \\\n')
            fout.write('-i {} \\\n'.format(input_file))
            fout.write('-p "cwl_lsf" \\\n')
            fout.write('-k {} \\\n'.format(cluster_num))
            fout.write('-o {} \\\n'.format(output_dir))
            fout.write('--dims-km {} \\\n'.format(dim))
            fout.write('-j ./config_cwlexec.json\n')


def calcu_ARI(author):
    root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'

    true_label_file = '{}/datasets/SilverStd/{}/{}_true_label.txt'.format(root_dir, author, author)
    true_label = pd.read_csv(true_label_file, delimiter='\t', header=0)
    # print(true_label)

    out_dir = '{}/outputs/SilverStd/{}/MICA_MDS/complete'.format(root_dir, author)
    for dim in range(1, 51):
        output_dir = glob.glob('{}/dim_{}/mica-*'.format(out_dir, dim))
        cluster_mem_file = '{}/cwl_lsf_k5_tsne_ClusterMem.txt'.format(output_dir[0])
        predict_label = pd.read_csv(cluster_mem_file, delimiter='\t', index_col=0)
        # print(predict_label)
        merged = true_label.merge(predict_label, left_on='sample', right_index=True)
        # print(merged)
        ari = adjusted_rand_score(merged['type'], merged['label'])
        print('{}\t{}'.format(dim, ari))
        # break


if __name__ == "__main__":
    author = 'Chung'
    # run_MDS(author)
    calcu_ARI(author)

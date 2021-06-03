#!/usr/bin/env python3

import os


def create_GE_cmd(root_dir, level, author):
    return


def create_MDS_cmd(root_dir, level, author):
    root_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA'
    test_dir = '{}/tests/SilverStd/{}/MICA_MDS/cmd_files'.format(root_dir, author)
    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)
    input_file = '{}/datasets/SilverStd/{}/{}_MICA_input.txt'.format(root_dir, author, author)
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
            fout.write('-k 8 \\\n')
            fout.write('-o {} \\\n'.format(output_dir))
            fout.write('--dims-km {} \\\n'.format(dim))
            fout.write('-j ./config_cwlexec.json\n')

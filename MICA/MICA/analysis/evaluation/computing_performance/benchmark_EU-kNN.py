#!/usr/bin/env python3

from pathlib import Path


def grid_search_EU_kNN_parameters(exp_file, root_out_dir):
    num_neighbors = range(10, 101, 10)
    num_jobs = range(10, 31, 5)

    for num_n in num_neighbors:
        for num_j in num_jobs:
            nw = num_j + 2
            out_dir = '{}/nne_{}_nj_{}'.format(root_out_dir, num_n, num_j)

            header = '''#BSUB -P MICA_GE\n#BSUB -q compbio\n#BSUB -oo {}/MICA_GE.out -eo {}/MICA_GE.err\n#BSUB -n {}
#BSUB -R "span[hosts=1]"\n#BSUB -M 3000'''.format(out_dir, out_dir, nw)

            Path(out_dir).mkdir(parents=True, exist_ok=True)
            cmd = '{}\nmica ge -i {} -o {} -ar 4.0 -nw {} -nne {}'.format(header, exp_file,
                                                                          out_dir, num_j, num_n)
            with open('{}/MICA_GE.sh'.format(out_dir), 'w') as fout:
                fout.write(cmd)


def pbmc20k():
    exp_file = '/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/PBMC_20k/' \
               'PBMC_20k_MICA_input.h5ad'
    root_out_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/' \
                   'computing_performance/EU-kNN'
    grid_search_EU_kNN_parameters(exp_file, root_out_dir)


if __name__ == "__main__":
    pbmc20k()

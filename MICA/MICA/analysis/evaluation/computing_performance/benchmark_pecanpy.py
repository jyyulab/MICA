#!/usr/bin/env python3

from pathlib import Path
import numpy as np


def grid_search_pecanpy_parameters(exp_file, root_out_dir):
    # walk_len = range(10, 121, 10)
    # n_walks = range(10, 121, 10)
    # window_size = range(5, 101, 5)
    p_hyper = np.arange(0.2, 3.1, 0.2)
    q_hyper = np.arange(0.2, 3.1, 0.2)
    num_jobs = range(10, 31, 5)

    for p in p_hyper:
        p = round(p, 2)
        for q in q_hyper:
            q = round(q, 2)
            for nw in num_jobs:
                out_dir = '{}/phyper_{}_qhyper_{}_nw_{}'.format(root_out_dir, p, q, nw)

                header = '''#BSUB -P MICA_GE\n#BSUB -q compbio\n#BSUB -oo {}/MICA_GE.out -eo {}/MICA_GE.err\n#BSUB -n {}
#BSUB -R "span[hosts=1]"\n#BSUB -M 3000'''.format(out_dir, out_dir, 32)

                Path(out_dir).mkdir(parents=True, exist_ok=True)
                cmd = '{}\nmica ge -i {} -o {} -ar 4.0 -hp {} -hq {} -nw {}'.format(header, exp_file, out_dir,
                                                                                    p, q, nw)
                with open('{}/MICA_GE.sh'.format(out_dir), 'w') as fout:
                    fout.write(cmd)


def pbmc20k():
    exp_file = '/research/rgs01/project_space/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/datasets/SilverStd/PBMC_20k/' \
               'PBMC_20k_MICA_input.h5ad'
    root_out_dir = '/research/projects/yu3grp/scRNASeq/yu3grp/LiangDing/MICA/tests/SilverStd/PBMC_20k/' \
                   'computing_performance/pecanpy/hyperparameter'
    grid_search_pecanpy_parameters(exp_file, root_out_dir)


if __name__ == "__main__":
    pbmc20k()

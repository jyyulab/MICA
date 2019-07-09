#!/usr/bin/env python3

import sys
import utils

trans = sys.argv[3].upper()

# def mds(fig_num, in_mat_file, max_dim, out_dir, out_file_name, perplexity=30, print_plot="True")

if trans == "MDS":
    utils.mds(
            sys.argv[2],  # in_mat_file
            int(sys.argv[4]),  # max_dim
            sys.argv[1],  # out_file_name
            dist_method=sys.argv[5],
            )
elif trans == "LPL":
    utils.lpl(
             sys.argv[2],
             int(sys.argv[4]),
             sys.argv[1],
             dist_method=sys.argv[5],
             )
elif trans == "PCA":
    utils.pca(
             sys.argv[2],
             int(sys.argv[4]),
             sys.argv[1],
             dist_method=sys.argv[5],
             )
elif trans == "LPCA":
    utils.lpl(
             sys.argv[2],
             int(sys.argv[4]),
             sys.argv[1],
             dist_method=sys.argv[5],
             )

    utils.pca(
             sys.argv[2],
             int(sys.argv[4]),
             sys.argv[1],
             dist_method=sys.argv[5],
             )
else:
    print("Transformation method not supported!")

print("[PREP DONE] Method used for dimension reduce: " + trans)

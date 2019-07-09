#!/usr/bin/env python3

import sys
import utils

project_name = sys.argv[1]
method = sys.argv[2]
mat_files = sys.argv[3:]

utils.merge_dist_mats(mat_files, project_name, method)

dist_file = project_name + "_dist.h5"

if method == "mi":
    utils.norm_mi_mat(dist_file, project_name)

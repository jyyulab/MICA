#!/usr/bin/env python3
import sys
import pandas as pd
import utils

print(sys.argv)

mat = pd.HDFStore(sys.argv[1])
metrics = str(sys.argv[2]).lower()

params = mat["params"]
mat1 = mat[params.loc["key1", 0]]
mat2 = mat[params.loc["key2", 0]]
mat.close()

# calc_mi_mat(mat, mat1, mat2, bins, m, key, out_dir, out_file_name)
# mat1, mat2, bins, m, MI_indx, path[2], project_name + MI_indx
utils.calc_distance_mat(mat1, mat2, params, method=metrics)

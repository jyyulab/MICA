#!/usr/bin/env python3
import os
import sys
import utils

print(sys.argv)

input_file = sys.argv[1]
out_name = sys.argv[2]
slice_unit = sys.argv[3]

# frame = pd.read_csv(input_file, sep=sep_sym, index_col=0).iloc[:, 0:]
# read_file(file, sep, header, out_dir, out_file_name, index_col="0")
utils.read_file(input_file, out_name)
# input file, # name of output

# fun: patch_file(df_file, out_file_name)
# prepare h5 files (whole)
h5_tmp = out_name + ".h5.tmp"
utils.patch_file(h5_tmp, out_name)
os.remove(h5_tmp)

# fun: slice_file(df_file,out_file_name, slice_size="1000")
h5_whole = out_name + ".whole.h5"
utils.slice_file(h5_whole, out_name, slice_unit)

# fun: calc_prep(in_file, project_name, mie_out_dir)
h5_sliced = out_name + ".sliced.h5"  # update input files
utils.calc_prep(h5_sliced, out_name)

print("[INFO] --> [MIE-PREP] Finished.")

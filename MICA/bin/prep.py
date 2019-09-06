#!/usr/bin/env python3
"""
This script is responsible for pre-processing HDF5-format files for computation.
"""

import os
from MICA.bin import utils
import argparse


def main():
    """Handles arguments and calls the driver function"""
    head_description='Slice data, name parameters in input file, and generate data files with pairs of slices'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Input file')
    parser.add_argument('-o', '--output-file', metavar='STR', required=True, help='Output file name')
    parser.add_argument('-s', '--slice-unit', type=int, metavar='INT', required=True, help='Size of data slices')
    args = parser.parse_args()

    prep(args.input_file, args.output_file, args.slice_unit)


def prep(input_file, out_name, slice_unit):
    """Calls utility functions to preprocess input file

    Reads file into HDF5-format, adds parameters, slices data in file, generates several files with
    different combinations of slices

    Args:
        input_file (str): path to input text-file
        out_name   (str): common rootname of generated output files
        slice_unit (int): size of each slice of cell data in input text-file
    """

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

if __name__ == '__main__':
    main()

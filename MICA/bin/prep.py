#!/usr/bin/env python3
""" This script is responsible for pre-processing HDF5-format files for computation.
"""

import os
import argparse
import logging

from MICA.lib import utils
from MICA.lib import preprocessing


def main():
    """ Handles arguments and calls the driver function. """
    head_description = 'Slice input matrix into small matrices, naming parameters, ' \
                       'and generate data files with pairs of slices'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Input h5ad (with an anndata object) '
                                                                                 'or txt file')
    parser.add_argument('-o', '--output-file', metavar='STR', required=True, help='Output file name')
    parser.add_argument('-s', '--slice-unit', type=int, metavar='INT', required=True, default=1000,
                        help='Size of sliced small matrix')
    args = parser.parse_args()

    prep(args.input_file, args.output_file, args.slice_unit)


def prep(input_file, out_name, slice_unit):
    """ Preprocess input file to create sliced matrices.

    Reads file into HDF5-format, adds parameters, slices data in file, generates several files with
    different combinations of slices.

    Args:
        input_file (str): path to input text-file
        out_name   (str): common rootname of generated output files
        slice_unit (int): size of each slice of cell data in input text-file
    """
    logging.basicConfig(level=logging.INFO)
    # utils.read_file(input_file, out_name)
    preprocessing.read_write_mat(input_file, out_name)

    # prepare h5 files (whole)
    h5_tmp = out_name + ".h5.tmp"
    utils.patch_file(h5_tmp, out_name)
    os.remove(h5_tmp)

    h5_whole = out_name + ".whole.h5"
    utils.slice_file(h5_whole, out_name, slice_unit)

    h5_sliced = out_name + ".sliced.h5"  # update input files
    utils.calc_prep(h5_sliced, out_name)

    logging.info('MICA-prep step completed successfully.')


if __name__ == '__main__':
    main()

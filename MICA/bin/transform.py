#!/usr/bin/env python3
import sys
import argparse
from MICA.bin import utils

def main():
    """Handles arguments and calls the driver function"""
    head_description = '' 
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Merged matrix file')
    parser.add_argument('-o', '--proj-name', metavar='STR', required=True, help='Project/output file name')
    parser.add_argument('-t', '--transform', metavar='STR', required=True, help='Transform method')
    parser.add_argument('-d', '--max-dimension', type=int, metavar='INT', required=True, help='Maximum number of dimensions in reduction')
    parser.add_argument('-m', '--metric', metavar='STR', required=True, help='Metric used in calculating distance')
    args=parser.parse_args()

    transform(args.input_file, args.proj_name, args.transform.upper(), args.max_dimension, args.metric)


def transform(infile, proj_name, trans, max_dim, metric):
    """Applies transform function to a matrix
    
    Args:
        infile    (str): path to input file
        proj_name (str): name of output file, ideally same as project name used in previous steps
        trans     (str): method used for transformation
        max_dir   (int): maximum numer of dimensions used in reduction
        method    (str): method used for distance calculation in previous steps
    """

    if trans == "MDS":
        utils.mds(infile, max_dim, proj_name, dist_method=metric)
    elif trans == "LPL":
        utils.lpl(infile, max_dim, proj_name, dist_method=metric)
    elif trans == "PCA":
        utils.pca(infile, max_dim, proj_name, dist_method=metric)
    elif trans == "LPCA":
        utils.lpl(infile, max_dim, proj_name, dist_method=metric)
        utils.pca(infile, max_dim, proj_name, dist_method=metric)
    else:
        print("Transformation method not supported!")

    print("[PREP DONE] Method used for dimension reduce: " + trans)

# def mds(fig_num, in_mat_file, max_dim, out_dir, out_file_name, perplexity=30, print_plot="True")

if __name__ == '__main__':
    main()

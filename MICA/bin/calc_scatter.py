#!/usr/bin/env python3
import pandas as pd
import argparse
from MICA.bin import utils

"""02_calc_scatter.py

This script is responsible for calculating calculating different metrics, used to compare gene 
expression levels in between cells.
"""

def main():
    """Handles arguments and calls the driver function"""
    head_description = 'Calculate either one of these metrics in between two slices of data: mutual information, euclidean distance, pearson or spearman correlations.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Input file')
    parser.add_argument('-m', '--metric', metavar='STR', required=True, help='Metric used in calculation')
    args=parser.parse_args()

    calc_scatter(args.input_file, args.metric.lower())


def calc_scatter(input_file, metric):
    """Calls calc_distance_mat utility function and calculates a metric in between cells that is chosen by the user

    Args:
        input_file (str): path to input HDF5-format file
        metric     (str): metric for calculation (mutual info, euclidean dist, pearson or spearman correlations
    """

    mat = pd.HDFStore(input_file)
    metrics = metric.lower()

    params = mat["params"]
    mat1 = mat[params.loc["key1", 0]]
    mat2 = mat[params.loc["key2", 0]]
    mat.close()

    # calc_mi_mat(mat, mat1, mat2, bins, m, key, out_dir, out_file_name)
    # mat1, mat2, bins, m, MI_indx, path[2], project_name + MI_indx
    utils.calc_distance_mat(mat1, mat2, params, method=metrics)

if __name__ == '__main__':
    main()

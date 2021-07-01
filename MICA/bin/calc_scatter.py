#!/usr/bin/env python3
import pandas as pd
import argparse
from MICA.lib import distance
from MICA.lib import utils

"""
This script is responsible for calculating calculating different metrics, used to compare gene 
expression levels in between cells.
"""


def main():
    """Handles arguments and calls the driver function"""
    head_description = 'Calculate either one of these metrics in between two slices of data: mutual information, ' \
                       'euclidean distance, pearson or spearman correlations.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Input file')
    parser.add_argument('-m', '--metric', metavar='STR', required=True, help='Metric used in calculation')
    args = parser.parse_args()

    calc_scatter(args.input_file, args.metric.lower())


def calc_scatter(input_file, metric):
    """ Calls calc_distance_mat utility function and calculates a metric in between cells that is chosen by the user

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
    utils.calc_distance_mat(mat1, mat2, params, method=metrics)
    # distance.calc_dis_mat_paras(mat1, mat2, params)


if __name__ == '__main__':
    main()

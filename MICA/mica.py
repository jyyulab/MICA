#!/usr/bin/env python3
# MICA entry point for both GE and MDS versions

import os
import sys
import logging
import time
import argparse
import numpy as np
import pandas as pd
import pathlib

from MICA.lib import neighbor_graph as ng
from MICA.lib import preprocessing as pp
from MICA.lib import dimension_reduction as dr
from MICA.lib import clustering as cl
from MICA.lib import visualize as vs


def main():
    head_description = 'MICA - Mutual Information-based Clustering Analysis tool.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every sub parser
    parent_parser = argparse.ArgumentParser(description='Parser for common arguments', add_help=False)
    parent_parser.add_argument('-i', '--input-file', metavar='FILE', required=True,
                               help='Path to an input file (h5ad file or tab-delimited text file)')
    parent_parser.add_argument('-o', '--output-dir', metavar='DIR', required=True, help='Path to final '
                                                                                        'output directory')
    parent_parser.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                               help='Visualization method UMAP or t-SNE (default: UMAP)')
    parent_parser.add_argument('-md', '--min-dist', metavar='FLOAT', required=False, default=0.6, type=float,
                               help='min_dist parameter in UMAP, minimum distance of points in the embedded space '
                                    '(default: 0.6)')
    parent_parser.add_argument('-pp', '--perplexity', metavar='INT', default=30,
                               help='Perplexity parameter for tSNE visualization')

    subparsers = parser.add_subparsers(title='subcommands', help='versions', dest='version')
    subparsers.required = True

    # Create a sub parser for running auto mode
    subparser_auto = subparsers.add_parser('auto', parents=[parent_parser], help='automatic version')
    subparser_auto.add_argument('-cc', '--cell-count', metavar='INT', required=False, default=5000, type=int,
                                help='Run MDS version if less than cell_count; otherwise, run GE version')

    # Create a sub parser for running graph embedding version
    subparser_ge = subparsers.add_parser('ge', parents=[parent_parser], help='graph embedding version')
    subparser_ge.add_argument('-dd', '--dr-dim', metavar='INT', required=False, default=20, type=int,
                              help='Number of dimensions to reduce to (default: 20)')
    subparser_ge.add_argument('-ir', '--min-resolution', metavar='FLOAT', required=False, default=0.4, type=float,
                              help='Determines the minimum size of the communities (default: 0.4)')
    subparser_ge.add_argument('-ar', '--max-resolution', metavar='FLOAT', required=False, default=3.4, type=float,
                              help='Determines the maximum size of the communities (default: 3.4)')
    subparser_ge.add_argument('-ss', '--step-size', metavar='FLOAT', required=False, default=0.4, type=float,
                              help='Determines the step size to sweep resolution from min_resolution to max_resolution '
                                   '(default: 0.4)')
    subparser_ge.add_argument('-nw', '--num-workers', metavar='INT', required=False, default=10, type=int,
                              help='Number of works to run in parallel (default: 10)')

    # Create a sub parser for running MDS version
    subparser_mds = subparsers.add_parser('mds', parents=[parent_parser], help='MDS version')
    subparser_mds.add_argument('-pn', '--project-name', metavar='STR', required=True, type=str,
                               help='Project name/ID.')
    subparser_mds.add_argument('-nc', '--num-clusters', metavar='INT', nargs='+', required=True, type=int,
                               help='Number of cluster to be specified in kmeans')
    subparser_mds.add_argument('-bs', '--bootstrap', metavar='INT', default=10, type=int,
                               help='Maximum number of iterations per dimension (default: 10)')
    subparser_mds.add_argument('-df', '--dist-func', metavar='STR', default="mi", type=str,
                               help='Method for distance matrix calculation [mi | euclidean | spearman | pearson]'
                                    '(default:mi)')
    subparser_mds.add_argument('-dm', '--dr-method', metavar='STR', default='MDS',
                               help='Transformation method used for dimension reduction '
                                    '[MDS | PCA | LPL | LPCA] (default: MDS)')
    subparser_mds.add_argument('-sn', '--slice-size', metavar='INT', default=1000, type=int,
                               help='Number of cells in each MI sub-matrix (default: 1000)')
    subparser_mds.add_argument('-tn', '--thread-number', metavar='INT', default=10, type=int,
                               help='Number of poolings used for multiple kmeans iterations,'
                                    'usually equals to iterations_km (default: 10)')
    subparser_mds.add_argument('-dk', '--dims-km', metavar='INT', nargs='+', default=[12, 13, 14, 15, 16, 17, 18, 19],
                               type=int, required=False,
                               help='Dimensions used in k-mean clustering, array inputs are supported '
                                    '(default: [12, 13, 14, 15, 16, 17, 18, 19])')
    subparser_mds.add_argument('-dp', '--dims-plot', metavar='INT', default=19, type=int,
                               help='Number of dimensions used in visualization (default: 19)')

    platform_subparsers = subparser_mds.add_subparsers(title='subcommands', help='platforms', dest='action')
    platform_subparsers.required = True

    # Create a sub parser for running cwltool
    subsubparser_local = platform_subparsers.add_parser('local', parents=[parent_parser],
                                                        help='run cwltool in a local workstation')
    subsubparser_local.add_argument('-sr', '--serial-run', help='run cwltool in serial mode', action='store_false')

    # Create a sub parser for running cwlexec
    subsubparser_lsf = platform_subparsers.add_parser('lsf', parents=[parent_parser],
                                                      help='run cwlexec in a IBM LSF interface')
    subsubparser_lsf.add_argument('-cj', '--config-json', metavar='FILE', required=True,
                                  help='LSF-specific configuration file in JSON format to be '
                                       'used for workflow execution')
    subsubparser_lsf.add_argument('-rr', '--rerun', metavar='STR', type=str,
                                  help='Rerun an exited workflow with the given workflow ID.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
    main()

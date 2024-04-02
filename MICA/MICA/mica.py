#!/usr/bin/env python3
# MICA entry point for both GE and MDS versions

import sys
import logging
import time
import argparse

from MICA.lib import preprocessing as pp
from MICA import mica_ge
from MICA import mica_mds


def main():
    head_description = 'MICA - Mutual Information-based Clustering Analysis tool.'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every sub parser
    parent_parser = argparse.ArgumentParser(description='Parser for common arguments', add_help=False)
    common_required = parent_parser.add_argument_group('common required arguments')
    common_required.add_argument('-i', '--input-file', metavar='FILE', required=True,
                                 help='Path to an input file (h5ad file or tab-delimited text file)')
    common_required.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                                 help='Path to final output directory')

    common_optional = parent_parser.add_argument_group('Optional arguments')
    common_optional.add_argument('-pn', '--project-name', metavar='STR', required=False, type=str,
                                 help='Project name/ID.')
    common_optional.add_argument('-nc', '--num-clusters', metavar='INT', nargs='+', required=False, type=int,
                                 help='Number of clusters to be specified in kmeans')
    common_optional.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                                 help='Visualization method UMAP or t-SNE (default: UMAP)')
    common_optional.add_argument('-md', '--min-dist', metavar='FLOAT', required=False, default=0.6, type=float,
                                 help='min_dist parameter in UMAP, minimum distance of points in the embedded space '
                                      '(default: 0.6)')
    common_optional.add_argument('-pp', '--perplexity', metavar='INT', default=30,
                                 help='Perplexity parameter for tSNE visualization')

    subparsers = parser.add_subparsers(title='subcommands', help='versions', dest='version')
    subparsers.required = True

    # Create a sub parser for running auto mode
    subparser_auto = subparsers.add_parser('auto', parents=[parent_parser], help='automatic version')
    add_required_auto = subparser_auto.add_argument_group('additional optional arguments')
    add_required_auto.add_argument('-cc', '--cell-count', metavar='INT', required=False, default=5000, type=int,
                                   help='Run MDS version if less than cell_count; otherwise, run GE version '
                                        '(default: 5000)')

    # Create a sub parser for running graph embedding version
    subparser_ge = subparsers.add_parser('ge', parents=[parent_parser], help='graph embedding version')
    mica_ge.add_ge_arguments(subparser_ge)

    # Create a sub parser for running MDS version
    subparser_mds = subparsers.add_parser('mds', parents=[parent_parser], help='MDS version')
    mica_mds.add_mds_arguments(subparser_mds)

    if len(sys.argv) == 1 or len(sys.argv) == 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if args.version == 'auto':
        start = time.time()
        logging.info('Start auto mode. Getting the number of cells ...')
        adata = pp.read_preprocessed_mat(args.input_file)
        frame = adata.to_df()
        end = time.time()
        runtime = end - start
        num_cells = frame.shape[0]
        logging.info('Done. {} cells. Runtime: {} seconds'.format(num_cells, runtime))

        if num_cells >= args.cell_count:    # Run GE version
            logging.info('More than {} cells. Use GE version...'.format(args.cell_count))
            mica_ge.add_ge_arguments(subparser_auto)
            args = parser.parse_args()
            mica_ge.mica_ge(args)
        else:                               # Run MDS version
            logging.info('Less than {} cells. Use MDS version...'.format(num_cells))
            mica_mds.add_mds_arguments(subparser_auto)
            args = parser.parse_args()
            if args.project_name is None:
                sys.exit('ArgumentError: argument -pn/--project-name is required')
            if args.num_clusters is None:
                sys.exit('ArgumentError: argument -nc/--num-clusters is required')
            mica_mds.mica_mds(args)
    elif args.version == 'ge':
        logging.info('Start GE mode...')
        mica_ge.mica_ge(args)
    elif args.version == 'mds':
        if args.platform == 'lsf' and args.config_json is None:
            sys.exit('Error: --config-json must be specified for lsf platform.')
        logging.info('Start MDS mode...')
        mica_mds.mica_mds(args)
    else:
        sys.exit('Error - invalid running version')


if __name__ == "__main__":
    main()

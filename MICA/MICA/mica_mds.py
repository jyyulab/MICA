#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shlex
import logging
import pathlib


def main():
    head_description = '''MICA - Mutual Information-based Clustering Analysis tool. This version uses a 
    multidimensional scaling method for dimension reduction.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input-file', metavar='FILE', required=True,
                          help='Path to an input file (h5ad file or tab-delimited text file)')
    required.add_argument('-o', '--output-dir', metavar='DIR', required=True, help='Path to final output directory')
    required.add_argument('-pn', '--project-name', metavar='STR', required=True, type=str,
                          help='Project name/ID.')
    required.add_argument('-nc', '--num-clusters', metavar='INT', nargs='+', required=True, type=int,
                          help='Number of clusters to be specified in kmeans')
    add_mds_arguments(parser)

    visualization = parser.add_argument_group('visualization arguments')
    visualization.add_argument('-vm', '--visual-method', metavar='STR', required=False, default='UMAP', type=str,
                               help='Visualization method UMAP or t-SNE (default: UMAP)')
    visualization.add_argument('-md', '--min-dist', metavar='FLOAT', required=False, default=0.6, type=float,
                               help='min_dist parameter in UMAP, minimum distance of points in the embedded space '
                                    '(default: 0.6)')
    visualization.add_argument('-pp', '--perplexity', metavar='INT', default=30,
                               help='Perplexity parameter for tSNE visualization')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mica_mds(args)


def add_mds_arguments(parser):
    platform = parser.add_argument_group('platform arguments')
    platform.add_argument('-pl', '--platform', metavar='STR', required=False, type=str, choices=['local', 'lsf'],
                          default='local', help='Platform local (cwltool) or lsf (cwlexec) (default: local)')
    platform.add_argument('-sr', '--serial-run', help='run cwltool in serial mode', action='store_false')
    platform.add_argument('-cj', '--config-json', metavar='FILE', required=False,
                          help='[Required if use lsf platform] LSF-specific configuration file in JSON format '
                               'to be used for workflow execution')
    platform.add_argument('-rr', '--rerun', metavar='STR', type=str, required=False,
                          help='Rerun an exited workflow with the given cwlexec workflow ID.')

    additional = parser.add_argument_group('additional arguments')
    additional.add_argument('-bs', '--bootstrap', metavar='INT', default=10, type=int,
                            help='Maximum number of iterations per dimension (default: 10)')
    additional.add_argument('-df', '--dist-func', metavar='STR', default="mi", type=str,
                            help='Method for distance matrix calculation [mi | euclidean | spearman | pearson]'
                                 '(default:mi)')
    additional.add_argument('-dm', '--dr-method', metavar='STR', default='MDS', required=False,
                            help='Transformation method used for dimension reduction '
                                 '[MDS | PCA | LPL | LPCA] (default: MDS)')
    additional.add_argument('-sn', '--slice-size', metavar='INT', default=1000, type=int, required=False,
                            help='Number of cells in each MI sub-matrix (default: 1000)')
    additional.add_argument('-tn', '--thread-number', metavar='INT', default=10, type=int, required=False,
                            help='Number of poolings used for multiple kmeans iterations,'
                                 'usually equals to iterations_km (default: 10)')
    additional.add_argument('-dk', '--dims-km', metavar='INT', nargs='+', default=[19], type=int, required=False,
                            help='Dimensions used in k-mean clustering, array inputs are supported '
                                 '(default: [19])')   # Use [12, 13, 14, 15, 16, 17, 18, 19] for robustness but slower
    additional.add_argument('-dp', '--dims-plot', metavar='INT', default=19, type=int, required=False,
                            help='Number of dimensions used in visualization (default: 19)')

    graph = parser.add_argument_group('graph clustering arguments')
    graph.add_argument('-nn', '--num-neighbor', metavar='INT', required=False, type=int,
                       help='Number of neighbors of a cell for building k-nearest neighbor graph')
    graph.add_argument('-dg', '--dims-graph', metavar='INT', nargs='+', default=19, required=False,
                       help='Dimensions used in graph clustering, array inputs are supported (default: 19)')
    return parser


def mica_mds(args):
    if args.platform == 'lsf' and args.config_json is None:
        sys.exit('Error: --config-json must be specified for lsf platform.')

    # to make executable and config findable
    installed_path = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] += (os.pathsep + installed_path + '/bin')
    cwl_path = installed_path + '/cwl'

    if args.num_neighbor:   # Use graph-based clustering
        cwl_script = 'mica_g.cwl'
    else:
        cwl_script = 'mica.cwl'

    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    if args.platform == 'local':
        fp_yml = create_input_yml(args)
        if args.serial_run:
            cmd = 'cwltool --leave-tmpdir --outdir {} {}/{} {}'.format(args.output_dir, cwl_path, cwl_script,
                                                                       fp_yml.name)
        else:
            cmd = 'cwltool --leave-tmpdir --parallel --preserve-environment HDF5_USE_FILE_LOCKING --leave-tmpdir ' \
                  '--outdir {} {}/{} {}'.format(args.output_dir, cwl_path, cwl_script, fp_yml.name)
    elif args.platform == 'lsf':
        os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
        if args.rerun:
            cmd = 'cwlexec -r {} -pe PATH -pe HDF5_USE_FILE_LOCKING -c {}'.format(
                   args.rerun, args.config_json)
        else:
            fp_yml = create_input_yml(args)
            cmd = 'cwlexec -pe PATH -pe HDF5_USE_FILE_LOCKING -c {} --outdir {} {}/{} {}'.format(
                   args.config_json, args.output_dir, cwl_path, cwl_script, fp_yml.name)
    else:
        sys.exit('Error - invalid platform.')

    logging.info(cmd)
    run_shell_command_call(cmd)
    logging.info('All done.')


def create_input_yml(args):
    """ Create input yml file. """
    with open(pathlib.PurePath(args.output_dir).joinpath('mica.yml'), 'w') as fp_yml:
        logging.info(fp_yml.name)
        contents = 'infile:\n  class: File\n  path: {}\n' \
                   'project_name: {}\n' \
                   'k: {}\n' \
                   'num_neighbor: {}\n' \
                   'visualization: {}\n' \
                   'dim_reduction: {}\n' \
                   'iterations_km: {}\n' \
                   'dims_km: {}\n' \
                   'dims: {}\n' \
                   'dims_plot: {}\n' \
                   'perplexity: {}\n' \
                   'min_dist: {}\n' \
                   'slice_size: {}\n' \
                   'thread_number: {}\n' \
                   'dist_metrics: {}\n'.format(os.path.abspath(args.input_file), args.project_name, args.num_clusters,
                                               args.num_neighbor, args.visual_method, args.dr_method.upper(),
                                               args.bootstrap, args.dims_km, args.dims_graph, args.dims_plot,
                                               args.perplexity, args.min_dist, args.slice_size, args.thread_number,
                                               args.dist_func)
        logging.info(contents)
        fp_yml.write(contents)
        fp_yml.flush()
        fp_yml.seek(0)
    return fp_yml


def run_shell_command_call(cmd):
    """ Wrapper of subprocess.check_call to take a cmd string as input
    Args:
        cmd (str): command to run
    """
    cmd_to_exec = shlex.split(cmd)
    subprocess.check_call(cmd_to_exec)


if __name__ == "__main__":
    main()

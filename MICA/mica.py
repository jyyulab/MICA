#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shlex
import logging
import tempfile
import pathlib
import json


def main():
    head_description = '''MICA [scMINER] is a scalable tool to perform unsupervised scRNA-seq clustering analysis.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every sub parser
    parent_parser = argparse.ArgumentParser(description='Parser for common arguments', add_help=False)

    parent_parser.add_argument('-i', '--exp', metavar='FILE', required=True,
                               help='Path to an expression matrix file, only tab separated files are allowed')

    parent_parser.add_argument('-p', '--project-name', metavar='NAME', required=True, type=str,
                               help='Project name/ID.')

    parent_parser.add_argument('-k', '--clusters', metavar='INT', nargs='+', required=True, type=int,
                               help='Number of cluster to be specified in kmeans')

    parent_parser.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                               help='Path to final output directory')

    parent_parser.add_argument('-b', '--bootstrap', metavar='INT', default=10, type=int,
                               help='Maximum number of iterations per dimension (default: 10)')

    parent_parser.add_argument('-dr', '--dim_reduction', metavar='STRING', default='MDS',
                               help='Transformation method used for dimension reduction '
                                    '[MDS | PCA | LPL | LPCA] (default: MDS)')

    parent_parser.add_argument('-v', '--visualization', metavar='STRING', default="tsne",
                               help='Visualization options to plot final clustering results: umap and tsne')

    parent_parser.add_argument('-pp', '--perplexity', metavar='INT', default=30,
                               help='Perplexity parameter for tSNE visualization')

    parent_parser.add_argument('-d', '--min-dist', metavar='FLOAT', default=0.01,
                               help='Minimum distance parameter for UMAP visualization (optional)')

    parent_parser.add_argument('-sn', '--slice-size', metavar='INT', default=1000,
                               help='Number of cells in each MI sub-matrix (default: 1000)')

    parent_parser.add_argument('-t', '--thread-number', metavar='INT', default=10,
                               help='Number of pooling used for multiple kmeans iterations,'
                                    'usually equals to iterations_km (default: 10)')

    parent_parser.add_argument('--dims-km', metavar='INT', nargs='+', default=[19],
                               help='Dimensions used in clustering, array inputs are supported (default: 19)')

    parent_parser.add_argument('--dims-plot', metavar='INT', default=19,
                               help='Number of dimensions used in visualization (default: 19)')

    parent_parser.add_argument('-r', '--resource', default=[12000, 16000, 20000],
                               help="Memory assigned to step merge and norm, dimension reduce and clustering "
                                    "(default: 12000, 16000, 20000")

    subparsers = parser.add_subparsers(title='Subcommands', help='platforms', dest='subcommand')
    subparsers.required = True

    # Create a sub parser for running cwltool
    subparser_local = subparsers.add_parser('local', parents=[parent_parser], help='run cwltool in a local workstation')
    subparser_local.add_argument('-s', '--serial', help='run cwltool in serial mode', action='store_true')

    # Create a sub parser for running cwlexec
    subparser_lsf = subparsers.add_parser('lsf', parents=[parent_parser], help='run cwlexec in a IBM LSF interface')
    subparser_lsf.add_argument('-q', '--queue', metavar='STR', required=True, help='LSF queue to submit the workflow')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    # to make executable and config findable
    installed_path = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] += (os.pathsep + installed_path + '/bin')
    cwl_path = installed_path + '/cwl'

    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(pathlib.PurePath(tmpdirname).joinpath('mica.yml'), 'w') as fp_yml:
            logging.info(fp_yml.name)
            contents = 'infile:\n  class: File\n  path: {}\n' \
                       'project_name: {}\n' \
                       'k: {}\n' \
                       'visualization: {}\n' \
                       'dim_reduction: {}\n' \
                       'iterations_km: {}\n' \
                       'dims_km: {}\n' \
                       'dims_plot: {}\n' \
                       'perplexity: {}\n' \
                       'min_dist: {}\n' \
                       'slice_size: {}\n' \
                       'thread_number: {}\n'.format(os.path.abspath(args.exp), args.project_name, args.clusters,
                                                    args.visualization, args.dim_reduction, args.bootstrap,
                                                    args.dims_km, args.dims_plot, args.perplexity, args.min_dist,
                                                    args.slice_size, args.thread_number)

            logging.info(contents)
            fp_yml.write(contents)
            fp_yml.flush()
            fp_yml.seek(0)

            if args.subcommand == 'local':
                if args.serial:
                    cmd = 'cwltool --outdir {} {}/mica.cwl {}'.format(args.output_dir,
                                                                                  cwl_path,
                                                                                  fp_yml.name)

                else:
                    cmd = 'cwltool --parallel --debug --outdir {} {}/mica.cwl {}'.format(args.output_dir,
                                                                                             cwl_path,
                                                                                             fp_yml.name)
            elif args.subcommand == 'lsf':
                with open(pathlib.PurePath(tmpdirname).joinpath('config.json'), 'w') as fp_config:
                    config_dict = {"queue": args.queue,
                                   "rerunnable": True,
                                   "steps": {
                                       "mergeAndnorm": {"res_req": "rusage[mem=" + str(args.resource[0]) + "]"},
                                       "dimension_reduce": {"res_req": "rusage[mem=" + str(args.resource[1]) + "]"},
                                       "clustering": {"res_req": "rusage[mem=" + str(args.resource[2]) + "]"}
                                   }
                                   }
                    logging.info(config_dict)
                    json.dump(config_dict, fp_config)
                    fp_config.flush()
                    fp_config.seek(0)
                    cmd = 'cwlexec -pe PATH -c {} --outdir {} ./MICA/cwl/mica.cwl {}'.format(
                        fp_config.name, args.output_dir, fp_yml.name)
            else:
                sys.exit('Error - invalid sub command.')
            logging.info(cmd)
            run_shell_command_call(cmd)

    logging.info('All done.')


def run_shell_command_call(cmd):
    """ Wrapper of subprocess.check_call to take a cmd string as input
    Args:
        cmd (str): command to run
    """
    cmd_to_exec = shlex.split(cmd)
    subprocess.check_call(cmd_to_exec)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
""" This module contains "auxiliary" functions.
"""

import os
import sys
import pathlib
import glob
import subprocess
import time
import datetime
import logging
import shlex


def run_shell_command(command_string):
    """ Executes a command and returns stdout, stderr, return_code.
        Input:
            - command_string: Command to be executed
        Output:
            - stdout: stdout of command as a single string.
            - stderr: stderr of command as a single string.
            - return_code: integer return code of command.
    """
    command = shlex.split(command_string)
    proc = subprocess.run(command, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    stdout = proc.stdout.decode('utf-8')
    stderr = proc.stderr.decode('utf-8')

    return_code = proc.returncode

    return stdout, stderr, return_code


def run_shell_command_call(cmd):
    """ Wrapper of subprocess.check_call to take a cmd string as input
    Args:
        cmd (str): command to run
    """
    cmd_to_exec = shlex.split(cmd)
    subprocess.check_call(cmd_to_exec)


def run_shell_command_output(cmd):
    """ Wrapper of subprocess.check_output to take a cmd string as input
    Args:
        cmd (str): command to run
    Returns:
        output of check_output
    """
    cmd_to_exec = shlex.split(cmd)
    return subprocess.check_output(cmd_to_exec).decode('utf-8')


def timeit(method):
    """ Decorator to calculate elapse time.
    """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print('%r %2.2f sec' % (method.__name__, te-ts), file=sys.stderr)
        return result
    return timed


def create_logger(log_root_dir, log_type, logger_name, **kwargs):
    """ Create a logger with a log file handler and a console handler.
    Args:
        log_root_dir (str): path to the root dir to save log files
        log_type (str): $SJCB_ENV_ROOT/logs/<type>
        logger_name (str): logger name
        optional keyword arguments:
            subtype (str): $SJCB_ENV_ROOT/logs/<type>/<subtype>
    Returns:
        A logger obj
    """
    subtype = kwargs.pop('subtype', None)
    year_str = str(datetime.datetime.now().strftime("%Y"))
    month_str = str(datetime.datetime.now().strftime("%m"))
    day_str = str(datetime.datetime.now().strftime("%d"))
    time_str = str(datetime.datetime.now().strftime("%H-%M-%S"))
    if not os.path.isdir(log_root_dir):
        sys.exit('Error - log root directory not found: {}'.format(log_root_dir))
    if subtype:
        log_dir = '{}/logs/{}/{}/{}/{}/{}'.format(log_root_dir, log_type, subtype, year_str, month_str, day_str)
    else:
        log_dir = '{}/logs/{}/{}/{}/{}'.format(log_root_dir, log_type, year_str, month_str, day_str)
    pathlib.Path(log_dir).mkdir(parents=True, exist_ok=True)
    log_file = '{}/log_{}.txt'.format(log_dir, time_str)
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_file)  # File handler
    fh.setLevel(logging.INFO)
    fh_formatter = logging.Formatter('%(asctime)s %(name)s:%(funcName)s:%(levelname)s:%(message)s')
    fh.setFormatter(fh_formatter)
    sh = logging.StreamHandler()  # Console handler
    sh.setLevel(logging.WARNING)
    sh_formatter = logging.Formatter('%(levelname)s:%(message)s')
    sh.setFormatter(sh_formatter)
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def read_config_sps(config_file):
    """ Read tab delimited config file with a sample, project, subproject per row and
        save to a list.
        Args:
            config_file (str): path to a config file
        Returns:
            list: a list of (sample, project, subproject)
    """
    sps_lst = []
    with open(config_file, 'r') as fin:
        while True:
            line = fin.readline().rstrip()
            if not line:
                break
            sps_lst.append(tuple(line.split('\t')))
    return sps_lst


def send_email(to_addrs, subj, body, from_addr):
    """ Send an email.
        Args:
            to_addrs (str): comma separated email addresses
            subj (str): subject
            body (str): email body
            from_addr (str): from email address
        Returns:
            None
    """
    cmd_str = 'send_mail.sh "{}" "{}" "{}" -r "{}"'.format(to_addrs, subj, body, from_addr)
    run_shell_command_call(cmd_str)


# Module wrapper and subcommand methods
def get_subcmds():
    """ Get all script names in the dir of the wrapper script and
        parse out the file names (as subcommands) and file types.
    Args:
        None
    Returns:
        (subcmd_list, subcmd:file_type dict)
    """
    pathname = os.path.dirname(sys.argv[0])
    abs_path = os.path.abspath(pathname)
    subcmd_paths = glob.glob(abs_path + '/*')
    subcmd_paths.remove(abs_path + '/' + os.path.basename(sys.argv[0]))
    subcmd_paths.sort()
    filenames = [os.path.basename(path) for path in subcmd_paths]
    filetype_dic = dict()
    subcmds = []
    for name in filenames:
        if name.find('.') != -1:
            tokens = name.split('.')
            subcmds.append(tokens[0])
            filetype_dic[tokens[0]] = tokens[1]
        else:
            subcmds.append(name)
            filetype_dic[name] = ''
    return subcmds, filetype_dic


def print_subcmds():
    """ Print all subcommands (file names without file types).
    """
    for subcmd in get_subcmds()[0]:
        print(subcmd)
    sys.exit(0)
    return 1


def subcmd_handler(file_type):
    """ Subcommand delegator.
    """
    if file_type:
        argLst = [sys.argv[1] + '.' + file_type]
    else:
        argLst = [sys.argv[1]]
    for arg in sys.argv[2:]:
        argLst.append(arg)
    subprocess.check_call(argLst)
    return 1

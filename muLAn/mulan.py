# -*-coding:Utf-8 -*
# =============================================================================
# MIT License
# 
# Copyright (c) 2014-2017 The muLAn project team
# Copyright (c) 2014-2017 Clément Ranc & Arnaud Cassan
# 
# See the LICENCE file
# =============================================================================

# --------------------------------------------------------------------
# Packages
# --------------------------------------------------------------------
import os
import sys
import importlib
import muLAn
# --------------------------------------------------------------------
# Fonctions
# --------------------------------------------------------------------
def add_topythonpath(path):
    if len(path)>0:
        sys.path.insert(0, path)
# --------------------------------------------------------------------
def check_slash(path):
    if len(path) > 0:
        if path[-1] != '/':
            path = '{:s}/'.format(path)
    return path
# --------------------------------------------------------------------
def gettxt_inifile(fn, kw, sep=':'):
    try:
        file = open(fn, 'r')
        result = [line.strip().split(sep)[1].strip() for line in
                file if line.strip().split(sep)[0].strip() == kw][-1]
        file.close()
    except IOError:
        sys.exit('File {:s} not found.'.format(fn))

    return result
# --------------------------------------------------------------------
def getpath_event():
    """Return the path of the event directory.

    :return: path of event directory
    :rtype: string
    """
    path = os.path.realpath(__file__)
    return '{:s}/'.format('/'.join(path.split('/')[:-1]))
# --------------------------------------------------------------------
def getpath_mulan(path_event):
    """Return root path of muLAn program.

    :param path_event: path of the event
    :type path_event: string
    :return: root path of muLAn
    :rtype: string
    """
    fn_as = '{:s}advancedsetup.ini'.format(path_event)
    path = gettxt_inifile(fn_as, 'Code')
    path = check_slash(path)
    return path
# --------------------------------------------------------------------
def getpath_localpackages(path_event):
    fn_as = '{:s}advancedsetup.ini'.format(path_event)
    path = gettxt_inifile(fn_as, 'PathLocalPythonPackages')
    path = check_slash(path)
    return path
# --------------------------------------------------------------------
def getint_verbose(path_event):
    fn_as = '{:s}setup.ini'.format(path_event)
    v = int(gettxt_inifile(fn_as, 'Verbose'))
    return v
# --------------------------------------------------------------------
def print_welcome(path, verbose=1):
    if verbose > 0:
        fn_logo = path + 'Logo.txt'
        try:
            file = open(fn_logo, 'r')
            logo = ''.join([line for line in file])
            file.close()
        except IOError:
            sys.exit('File {:s} not found.'.format(fn_logo))

        legend = '\033[90mMicro-Lensing Analysis software\033[0m\n'
        version = '\033[90mVersion '+muLAn.__version__

        welcome = '\033[1m\033[34m{:s}\033[0m{:s}{:s}\033[0m'.format(logo, legend, version)
        print welcome
# --------------------------------------------------------------------
def bash_command(cmd):
    proc = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    proc.wait()
# --------------------------------------------------------------------
def check_packages(path_mulan, verbose=0):
    try:
        from distutils.version import LooseVersion
    except ImportError:
        txt = '\033[1m\033[31mThe required package distutils is not installed.\033[0m'
        txt += '\nPlease install it or run the command'
        txt += '\n   $ pip install -r Requirements.txt'
        txt += '\nfrom the following directory'
        txt += '\n   ' + path_mulan
        txt += '\nto install all the required packages.'
        txt += '\n\033[1m\033[31mmuLAn stopped.\033[0m'
        sys.exit(txt)

    fn_requirements = '{:s}Requirements.txt'.format(path_mulan)
    try:
        file = open(fn_requirements, 'r')
        result = [line.replace('\n', '').replace(' ', '').strip()
                  for line in file if line.replace('\n', '') != '']
        file.close()
    except IOError:
        sys.exit('File {:s} not found.'.format(fn_requirements))

    packages = []
    versions = []
    for a in result:
        end = len([True for l in a if l.isalpha()])
        if end > 0:
            packages.append(a[:end])
            versions.append(a[end:].replace('=', '').replace('>', '').replace('<', ''))
            if versions[-1] == '': versions[-1] = 'None'

    for p, v in zip(packages, versions):
        try:
            if p=='GetDist':
                mod = importlib.import_module(p.lower())
            else:
                mod = importlib.import_module(p)
        except ImportError:
            txt = '\033[1m\033[31mThe required package {:s} is not installed.\033[0m'.format(p)
            txt += '\nPlease install it or run the command'
            txt += '\n   $ pip install -r Requirements.txt'
            txt += '\nfrom the following directory'
            txt += '\n   ' + path_mulan
            txt += '\nto install all the required packages.'
            txt += '\n\033[1m\033[31mmuLAn stopped.\033[0m'
            sys.exit(txt)

        try:
            v_computer = mod.__version__
        except:
            v_computer = 'None'

        if (v_computer != 'None') & (v != 'None'):
            if LooseVersion(v) > LooseVersion(v_computer):
                txt = 'An earlier version of the package {:s} is required.'.format(p)
                txt += '\n   Installed: {:s}   Required: {:s}'.format(v_computer, v)
                txt += '\n\033[1m\033[31mmuLAn stopped.\033[0m'
                sys.exit(txt)
            if (LooseVersion(v) < LooseVersion(v_computer)) & verbose > 0:
                txt = 'An earlier version of the package {:s} is installed.'.format(p)
                txt += '\n   Installed: {:s}   Required: {:s}'.format(v_computer, v)
                print txt
# ----------------------------------------------------------------------
def run_mulan(path_event, options):
    mul = importlib.import_module('sequential')
    mul.run_sequence(path_event, options)
# ----------------------------------------------------------------------
def sort(path_event):
    # Sort exploration results
    order_ChainsResults = importlib.import_module('order_ChainsResults')
    order_ChainsResults.order(path_event)
# ----------------------------------------------------------------------
from muLAn.packages.order_ChainsResults import order
def sortmodels():
    ici = os.getcwd()
    path_event = os.path.realpath(ici)+'/'
    order(path_event)
# ----------------------------------------------------------------------
def stop():
    # Safe Emergency Stop
    ici = os.getcwd()
    path_event = os.path.realpath(ici)+'/'
    if os.path.exists(path_event + '.emergencystop'):
        os.remove(path_event + '.emergencystop')
        file = open(path_event + '.emergencystop', 'w')
        file.write('1')
        file.close()
    fn_lock = '{:s}.lock'.format(path_event) # Remove .lock file
    if os.path.exists(fn_lock):
        os.remove(fn_lock)
# ----------------------------------------------------------------------
def run():
    # Event (local) path
    ici = os.getcwd()
    path_event = os.path.realpath(ici)+'/'
    verbose = getint_verbose(path_event)
    path_mulan = getpath_mulan(path_event)
    print_welcome(path_mulan, verbose=verbose)
    # Add local packages
    path_localpackages = getpath_localpackages(path_event)
    add_topythonpath(path_localpackages)
    module_path = '{:s}packages'.format(path_mulan)
    add_topythonpath(module_path)
    # Check packages
    ##removed## check_packages(path_mulan, verbose=verbose)

    # Load standard packages (after adding path to local packages)
    argparse = importlib.import_module('argparse')
    subprocess = importlib.import_module('subprocess')
    np = importlib.import_module('numpy')

    # Command line options - TO BE REMOVED IN FORTHCOMING VERSION
    text = 'The command line options below overwrite the configuration file.'
    parser = argparse.ArgumentParser(prog='python mulan.py', description=text)
    parser.add_argument('-a', '--archive', nargs=1, type=str, default=['None'], help='Replace <ARCHIVE> by the name of the archive.')
    parser.add_argument('-f', '--fit', action='store_true', help='Ask muLAn to fit the data.')
    parser.add_argument('--nchains', nargs=1, type=int, default=[-1], help='Replace <NCHAINS> by the number of chains.')
    parser.add_argument('--ncores', nargs=1, type=int, default=[-1], help='Replace <NCORES> by the number of cores.')
    parser.add_argument('-o', '--optimize', action='store_true', help='Optimize inputs/outputs using pickle package when possible.')
    parser.add_argument('-p', '--plot', nargs=1, type=str, default=['None'], help='Ask muLAn to plot the model '
            + 'number <PLOT> (default: best fitting model). For several models, use coma separator without space '
            + '(e.g. -p 1,5,7) or a dash without space for a range of models (e.g. -p 1-5).')
    parser.add_argument('--resume', action='store_true', help='Resume an exploration.')
    parser.add_argument('-s', '--sort', action='store_true', help='Sort the samples and create a table.')
    parser.add_argument('-sn', '--sortno', action='store_true', help='Skip sort stage when running muLAn.')
    parser.add_argument('--stop', action='store_true', help='Stop the exploration and save properly the results.')
    parser.add_argument('-v', '--verbose', nargs=1, type=int, choices=range(0, 6), default=[-1], help='Choose a verbose level.')
    opts = parser.parse_args()
    options = dict()
    [options.update({key: getattr(opts, key)}) for key in vars(opts)]
    [options.update({key: np.atleast_1d(options[key])[0]}) for key in options if len(np.atleast_1d(options[key]))==1]
    if options['verbose'] > -1:
        verbose = options['verbose']

    # Add muLAn packages and modules
    sys.path.insert(0, path_mulan)


    cond = not options['sort'] and not options['stop']
    if cond:
        run_mulan(path_event, options)






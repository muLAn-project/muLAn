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
import muLAn.packages.sequential as sequential 
# --------------------------------------------------------------------
# Fonctions
# --------------------------------------------------------------------
def getlogo():
    logo =\
            "                 _       ___        \n"\
          + "                | |     / _ \       \n"\
          + " _ __ ___  _   _| |    / /_\ \_ __  \n"\
          + "| '_ ` _ \| | | | |    |  _  | '_ \ \n"\
          + "| | | | | | |_| | |____| | | | | | |\n"\
          + "|_| |_| |_|\__,_\_____/\_| |_/_| |_|\n"
    return logo
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
def print_welcome(verbose=1):
        logo = getlogo()
        legend = '\033[90mMicro-Lensing Analysis software\033[0m\n'
        version = '\033[90mVersion {:s}'.format(muLAn.__version__)

        welcome = '\033[1m\033[34m{:s}\033[0m{:s}{:s}\033[0m'.format(logo, legend, version)
        print welcome
# --------------------------------------------------------------------
def bash_command(cmd):
    proc = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    proc.wait()
# --------------------------------------------------------------------
def run_mulan(path_event, options):
    # mul = importlib.import_module('sequential')
    # mul.run_sequence(path_event, options)
    sequential.run_sequence(path_event, options)
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
def run(options=dict()):
    # Event (local) path
    ici = os.getcwd()
    path_event = os.path.realpath(ici)+'/'
    verbose = getint_verbose(path_event)
#    path_mulan = getpath_mulan(path_event)
    print_welcome(verbose=verbose)

    # Load standard packages (after adding path to local packages)
#    argparse = importlib.import_module('argparse')
#    subprocess = importlib.import_module('subprocess')
#    np = importlib.import_module('numpy')

    # Add muLAn packages and modules
#    sys.path.insert(0, path_mulan)

#    cond = not options['sort'] and not options['stop']
#    if cond:
    run_mulan(path_event, options)






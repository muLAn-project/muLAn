# -*-coding:Utf-8 -*
# ====================================================================
# MIT License
# 
# Copyright (c) 2014-2017 The muLAn project team
# Copyright (c) 2014-2017 Clément Ranc & Arnaud Cassan
# 
# See the LICENCE file
# ====================================================================

# --------------------------------------------------------------------
# Packages
# --------------------------------------------------------------------
import os
import sys
import importlib
import muLAn
import muLAn.packages.sortmodels as mulansort
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
    ici = os.getcwd()
    path_event = check_slash(os.path.realpath(ici))
    return path_event
# --------------------------------------------------------------------
def getint_verbose(path_event):
    fn_as = '{:s}setup.ini'.format(path_event)
    v = int(gettxt_inifile(fn_as, 'Verbose'))
    return v
# --------------------------------------------------------------------
def print_welcome(verbose=1):
        logo = getlogo()
        # legend = '\033[90mMicro-Lensing Analysis software\033[0m\n'
        legend = '\033[36mgravitational MICROlensing Analysis code\033[0m\n'
        version = '\033[36mVersion {:s}'.format(muLAn.__version__)
        welcome = '\033[1m\033[34m{:s}\033[0m{:s}{:s}\033[0m'.format(logo, legend, version)
        print(welcome)
# ----------------------------------------------------------------------
def check_options(options):
    if not 'sortno' in options:
        options.update({'sortno': False})
    if not 'plot' in options:
        options.update({'plot': None})
    if not 'fit' in options:
        options.update({'fit': None})
    if not 'archive' in options:
        options.update({'archive': None})
    if not 'ncores' in options:
        options.update({'ncores': None})
    if not 'nchains' in options:
        options.update({'nchains': None})
    if not 'resume' in options:
        options.update({'resume': None})
    if not 'verbose' in options:
        options.update({'verbose': None})
    if not 'optimize' in options:
        options.update({'optimize': None})
# --------------------------------------------------------------------
def run_mulan(path_event, options):
    sequential.run_sequence(path_event, options)
# ----------------------------------------------------------------------
def sortmodels():
    mulansort.McmcFiles().sort()
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
def run(options=None):

    # Welcoming message depending on verbose option 
    ici = os.getcwd()
    path_event = check_slash(os.path.realpath(ici))
    verbose = getint_verbose(path_event)
    print_welcome(verbose=verbose)

    # Use a default options dictionary if not provided by user
    if options == None: options = dict()
    check_options(options)

    # Run the modeling code
    run_mulan(path_event, options)


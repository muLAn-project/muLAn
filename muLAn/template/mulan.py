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
import argparse
import subprocess
import numpy as np
import muLAn.mulan as mulan



if __name__ == '__main__':

    # Command line options
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

    if options['stop']:
        mulan.stop()

    if options['sort']:
        mulan.sortmodels()

    if (not options['sort'] and not options['stop']):
        mulan.run(options=options)






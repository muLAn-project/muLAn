# -*-coding:Utf-8 -*
"""muLAn (gravitational MICROlensing Analysis code) run commands"""

# Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan
# Distributed under the terms of the MIT license
#
# This module is part of software:
#       muLAn: gravitational MICROlensing Analysis code
#       https://github.com/muLAn-project/muLAn

import argparse
import numpy as np
import muLAn.mulan as mulan
from muLAn.utils import muLAnFig
from muLAn.utils.muLAnCleanData import cleandata

# --------------------------------------------------------------------
# USER: CHOOSE RUN MODE
# --------------------------------------------------------------------
#   ComLine: boolean
#       True:  Execute muLAn from command line
#       False: Use 'setupe.ini' options and runs listed commands.
#              This mode further allows to execute utility routines such
#              as data cleaning of creating output pdf figures.
ComLine = False

# --------------------------------------------------------------------
# Execute muLAn with following instructions and 'setup.ini' options
# --------------------------------------------------------------------
if not ComLine:

    ### clean data?
#    cleandata('Data/data.dat')

    ### run muLAn using 'setup.ini' options?
    mulan.run()

    ### create a figure for publication?
#    help(muLAnFig)
#    fig = muLAnFig.figure()
#    fig.plot()
#    fig.addinset_caustics([0.2, 0.7, 0.2, 0.2])
#    fig.addinset_lightcurve([0.2, 0.4, 0.2, 0.2])
#    fig.save('Plots/figure.pdf')
#    fig.show()

# --------------------------------------------------------------------
# Execute muLAn from command line
# --------------------------------------------------------------------
else:
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

    # Actions
    if options['stop']:
        mulan.stop()

    if options['sort']:
        mulan.sortmodels()

    if (not options['sort'] and not options['stop']):
        mulan.run(options=options)


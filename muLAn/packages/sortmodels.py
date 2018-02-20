# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   External libraries
# ----------------------------------------------------------------------
import sys
import os
# ----------------------------------------------------------------------
#   Packages
# ----------------------------------------------------------------------
import glob
import copy
import pandas as pd
import muLAn
import muLAn.packages.general_tools as gtools
import numpy as np
import ConfigParser as cp
# ----------------------------------------------------------------------
#   CLASS
# ----------------------------------------------------------------------
class printoption:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'

    reset = '\033[0m'
    bright = '\033[1m'
    dim = '\033[2m'
    underscore = '\033[4m'
    blink = '\033[5m'
    reverse = '\033[7m'
    hidden = '\033[8m'

    level0 = "\033[1m\033[31m"
    level1 = "\033[1m"
    good = "\033[32m"


# ====================================================================
#  Class
# ====================================================================
class McmcFiles:
    """A class to read MCMC output files and sort results.
    
    :param str path: Path to the MCMC files.
    """
    def __init__(self, path=None):
        self._path = path

        if path==None:
            self._getconfig()

    def _getconfig(self):
        """Load the configuration files *.ini."""

        # Path of the event
        path_event = muLAn.mulan.getpath_event()

        # Configuration files
        fname_setup = "{:s}setup.ini".format(path_event)
        fname_advanced = "{:s}advancedsetup.ini".format(path_event)

        # Load configuration files
        cfgsetup = cp.SafeConfigParser()
        cfgsetup.read([fname_setup, fname_advanced])
        cfgobs = cp.SafeConfigParser()
        cfgobs.read(path_event + 'observatories.ini')

        # Add the path to the configuration
        cfgsetup.set('FullPaths', 'Event', path_event)

        # Check the paths
        if cfgsetup.get('FullPaths', 'Code').replace(" ", "") != "":
            if cfgsetup.get('FullPaths', 'Code')[-1] != '/':
                cfgsetup.set('FullPaths', 'Code', cfgsetup.get('FullPaths', 'Code') + '/')
            if cfgsetup.get('FullPaths', 'Event')[-1] != '/':
                cfgsetup.set('FullPaths', 'Event', cfgsetup.get('FullPaths', 'Event') + '/')
        if cfgsetup.get('RelativePaths', 'Data')[-1] != '/':
            cfgsetup.set('RelativePaths', 'Data', cfgsetup.get('RelativePaths', 'Data') + '/')
        if cfgsetup.get('RelativePaths', 'Plots')[-1] != '/':
            cfgsetup.set('RelativePaths', 'Plots', cfgsetup.get('RelativePaths', 'Plots') + '/')
        if cfgsetup.get('RelativePaths', 'Chains')[-1] != '/':
            cfgsetup.set('RelativePaths', 'Chains', cfgsetup.get('RelativePaths', 'Chains') + '/')
        if cfgsetup.get('RelativePaths', 'Outputs')[-1] != '/':
            cfgsetup.set('RelativePaths', 'Outputs', cfgsetup.get('RelativePaths', 'Outputs') + '/')
        if cfgsetup.get('RelativePaths', 'Archives')[-1] != '/':
            cfgsetup.set('RelativePaths', 'Archives', cfgsetup.get('RelativePaths', 'Archives') + '/')
        if cfgsetup.get('RelativePaths', 'ModelsHistory')[-1] != '/':
            cfgsetup.set('RelativePaths', 'ModelsHistory', cfgsetup.get('RelativePaths', 'ModelsHistory') + '/')

        self._cfgsetup = cfgsetup

    def sort(self):

        cfgsetup = self._cfgsetup
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')

        fnames_chains = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*-c*.txt")
        fnames_chains_exclude = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*g*.txt")

        temp =[]
        for a in fnames_chains:
            if (a in fnames_chains_exclude)==False:
                temp.append(a)

        fnames_chains = copy.deepcopy(temp)
        del temp, fnames_chains_exclude

        nb_chains = len(fnames_chains)

        if nb_chains!=0:

            samples_file = dict(
                {'chi2': [], 't0': [], 'u0': [], 'tE': [], 'rho': [], \
                 'gamma': [], 'piEE': [], 'piEN': [], 's': [], 'q': [], \
                 'alpha': [], 'dalpha': [], 'ds': [], 'chain': [], 'fullid': [],\
                 'date_save': [], 'time_save': [], 'id': [], 'accrate': [],\
                 'chi2/dof': []})

            # Test if an history already exist
            filename_history = cfgsetup.get('FullPaths', 'Event')\
                + cfgsetup.get('RelativePaths', 'ModelsHistory')\
                + cfgsetup.get('Controls', 'Archive')\
                + '-ModelsSummary.txt'

            if os.path.exists(filename_history):
            # if 0:
                file = open(filename_history, 'r')
                for line in file:
                    params_model = line

                    if params_model[0] == '#':
                        continue

                    try:
                        samples_file['fullid'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][0]))
                        samples_file['t0'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][1]))
                        samples_file['u0'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][2]))
                        samples_file['tE'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][3]))
                        samples_file['rho'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][4]))
                        samples_file['gamma'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][5]))
                        samples_file['piEN'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][6]))
                        samples_file['piEE'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][7]))
                        samples_file['s'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][8]))
                        samples_file['q'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][9]))
                        samples_file['alpha'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][10]))
                        samples_file['dalpha'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][11]))
                        samples_file['ds'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][12]))
                        samples_file['chi2'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][13]))
                        samples_file['chi2/dof'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][14]))
                        samples_file['accrate'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][15]))
                        samples_file['date_save'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][16]))
                        samples_file['time_save'].append(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][17])
                        samples_file['chain'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][18]))
                        samples_file['id'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][19]))
                    except:
                        text = "\n\033[1m\033[91mThe file\033[0m\n" + filename_history\
                               + "\n\033[1m\033[91mis corrupted. muLAn killed.\033[0m"
                        sys.exit(text)
                file.close()

            # Read on the chains
            if nb_chains > 0:
                for i in xrange(nb_chains):

                    file = open(fnames_chains[i], 'r')
                    for line in file:
                        params_model = line

                        if params_model[0] == '#':
                            continue

                        try:
                            samples_file['id'].append(int(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][0]))
                            samples_file['t0'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][1]))
                            samples_file['u0'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][2]))
                            samples_file['tE'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][3]))
                            samples_file['rho'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][4]))
                            samples_file['gamma'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][5]))
                            samples_file['piEN'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][6]))
                            samples_file['piEE'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][7]))
                            samples_file['s'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][8]))
                            samples_file['q'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][9]))
                            samples_file['alpha'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][10]))
                            samples_file['dalpha'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][11]))
                            samples_file['ds'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][12]))
                            samples_file['chi2'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][13]))
                            samples_file['accrate'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][14]))
                            samples_file['date_save'].append(int(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][15]))
                            samples_file['time_save'].append(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][16])
                            samples_file['chi2/dof'].append(float(
                                [a for a in (params_model.split('\n')[0].split(' ')) if
                                 (a != '')][17]))
                            samples_file['chain'].append(int(fnames_chains[i][-8:-4]))
                            samples_file['fullid'].append(-1)
                        except:
                            text = "\n\033[1m\033[91mThe file\033[0m\n" + "\033[1m\033[91m" + fnames_chains[i]\
                                   + "\033[0m\n\033[1m\033[91mis corrupted. muLAn killed.\033[0m"
                            sys.exit(text)
                    file.close()

            # Order the models
            chi2_min = np.min(samples_file['chi2'])
            samples_file.update(
                {'dchi2': samples_file['chi2'] - chi2_min})

            # Remove duplicates
            results = pd.DataFrame({})
            for key in samples_file:
                results[key] = samples_file[key]

            results_sorted = results.sort_values(['dchi2', 'fullid'], ascending=[1, 0]).drop_duplicates(
                    subset=['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds', 'chi2'])

            # Give a unique ID to models
            id_start = np.max(results_sorted['fullid']) + 1
            if id_start == 0 : id_start = 1
            cond = results_sorted['fullid'] == -1
            results_sorted.loc[cond, 'fullid'] = id_start + np.arange(cond.sum())

            # Save new file in csv with exponential
            filename_history = filename_history[:-3] + 'csv'
            file = open(filename_history, 'w')
            format = '#{:},'.format('UniqueID')\
                + '{:},'.format('t0')\
                + '{:},'.format('u0')\
                + '{:},'.format('tE')\
                + '{:},'.format('rho')\
                + '{:},'.format('gamma')\
                + '{:},'.format('piEN')\
                + '{:},'.format('piEE')\
                + '{:},'.format('s')\
                + '{:},'.format('q')\
                + '{:},'.format('alpha')\
                + '{:},'.format('dalpha')\
                + '{:},'.format('ds')\
                + '{:},'.format('chi2')\
                + '{:},'.format('chi2/dof')\
                + '{:},'.format('accrate')\
                + '{:}'.format('chain')\
                + '\n'
            file.write(format)

            for i in xrange(len(results_sorted)):
                format = '{:},'.format(results_sorted['fullid'].values[i])\
                    + '{:.10e},'.format(results_sorted['t0'].values[i])\
                    + '{:.10e},'.format(results_sorted['u0'].values[i])\
                    + '{:.10e},'.format(results_sorted['tE'].values[i])\
                    + '{:.10e},'.format(results_sorted['rho'].values[i])\
                    + '{:.10e},'.format(results_sorted['gamma'].values[i])\
                    + '{:.10e},'.format(results_sorted['piEN'].values[i])\
                    + '{:.10e},'.format(results_sorted['piEE'].values[i])\
                    + '{:.10e},'.format(results_sorted['s'].values[i])\
                    + '{:.10e},'.format(results_sorted['q'].values[i])\
                    + '{:.10e},'.format(results_sorted['alpha'].values[i])\
                    + '{:.10e},'.format(results_sorted['dalpha'].values[i])\
                    + '{:.10e},'.format(results_sorted['ds'].values[i])\
                    + '{:.10e},'.format(results_sorted['chi2'].values[i])\
                    + '{:.10e},'.format(results_sorted['chi2/dof'].values[i])\
                    + '{:.3f},'.format(results_sorted['accrate'].values[i])\
                    + '{:}'.format(results_sorted['chain'].values[i])\
                    + '\n'
                file.write(format)
            file.close()






# ----------------------------------------------------------------------
#   Functions
# ----------------------------------------------------------------------
def communicate(cfg, verbose, text, opts=False, prefix=False, newline=False, tab=False):
    if cfg.getint('Modelling', 'Verbose') >= verbose:
        if prefix:
            text = "[muLAn] " + text
        if opts!=False:
            text2=''
            for a in opts:
                text2 = text2 + a
                text = text2 + text + printoption.reset
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print text
        else:
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print text
# ----------------------------------------------------------------------
def unpack_options(cfgsetup, level0, level1):
    options = [a.strip() for a in cfgsetup.get(level0, level1).split(',')]
    del a, cfgsetup, level0, level1
    return options
# ----------------------------------------------------------------------
def order(full_path):
    # ----------------------------------------------------------------------
    #   Configuration files
    # ----------------------------------------------------------------------
    # Load configuration files
    cfgsetup = cp.SafeConfigParser()
    cfgsetup.read([full_path + 'setup.ini', full_path + 'advancedsetup.ini'])

    cfgobs = cp.SafeConfigParser()
    cfgobs.read(full_path + 'observatories.ini')

    # Add the path to the configuration
    cfgsetup.set('FullPaths', 'Event', full_path)

    # Check the paths
    if cfgsetup.get('FullPaths', 'Code')[-1]!='/':
        cfgsetup.set('FullPaths', 'Code', cfgsetup.get('FullPaths', 'Code') + '/')
    if cfgsetup.get('FullPaths', 'Event')[-1]!='/':
        cfgsetup.set('FullPaths', 'Event', cfgsetup.get('FullPaths', 'Event') + '/')
    if cfgsetup.get('RelativePaths', 'Data')[-1]!='/':
        cfgsetup.set('RelativePaths', 'Data', cfgsetup.get('RelativePaths', 'Data') + '/')
    if cfgsetup.get('RelativePaths', 'Plots')[-1]!='/':
        cfgsetup.set('RelativePaths', 'Plots', cfgsetup.get('RelativePaths', 'Plots') + '/')
    if cfgsetup.get('RelativePaths', 'Chains')[-1]!='/':
        cfgsetup.set('RelativePaths', 'Chains', cfgsetup.get('RelativePaths', 'Chains') + '/')
    if cfgsetup.get('RelativePaths', 'Archives')[-1]!='/':
        cfgsetup.set('RelativePaths', 'Archives', cfgsetup.get('RelativePaths', 'Archives') + '/')
    if cfgsetup.get('RelativePaths', 'ModelsHistory')[-1]!='/':
        cfgsetup.set('RelativePaths', 'ModelsHistory', cfgsetup.get('RelativePaths', 'ModelsHistory') + '/')

    # Order the models

    text = "Sort parameters"
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    #   Initialisation
    # ------------------------------------------------------------------

    path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')

    fnames_chains = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*-c*.txt")
    fnames_chains_exclude = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*g*.txt")

    temp =[]
    for a in fnames_chains:
        if (a in fnames_chains_exclude)==False:
            temp.append(a)

    fnames_chains = copy.deepcopy(temp)
    del temp, fnames_chains_exclude

    nb_chains = len(fnames_chains)

    if nb_chains!=0:

        samples_file = dict(
            {'chi2': [], 't0': [], 'u0': [], 'tE': [], 'rho': [], \
             'gamma': [], 'piEE': [], 'piEN': [], 's': [], 'q': [], \
             'alpha': [], 'dalpha': [], 'ds': [], 'chain': [], 'fullid': [],\
             'date_save': [], 'time_save': [], 'id': [], 'accrate': [],\
             'chi2/dof': []})

        # Test if an history already exist
        filename_history = cfgsetup.get('FullPaths', 'Event')\
            + cfgsetup.get('RelativePaths', 'ModelsHistory')\
            + cfgsetup.get('Controls', 'Archive')\
            + '-ModelsSummary.txt'

        if os.path.exists(filename_history):
        # if 0:
            file = open(filename_history, 'r')
            for line in file:
                params_model = line

                if params_model[0] == '#':
                    continue

                try:
                    samples_file['fullid'].append(int(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][0]))
                    samples_file['t0'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][1]))
                    samples_file['u0'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][2]))
                    samples_file['tE'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][3]))
                    samples_file['rho'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][4]))
                    samples_file['gamma'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][5]))
                    samples_file['piEN'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][6]))
                    samples_file['piEE'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][7]))
                    samples_file['s'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][8]))
                    samples_file['q'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][9]))
                    samples_file['alpha'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][10]))
                    samples_file['dalpha'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][11]))
                    samples_file['ds'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][12]))
                    samples_file['chi2'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][13]))
                    samples_file['chi2/dof'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][14]))
                    samples_file['accrate'].append(float(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][15]))
                    samples_file['date_save'].append(int(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][16]))
                    samples_file['time_save'].append(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][17])
                    samples_file['chain'].append(int(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][18]))
                    samples_file['id'].append(int(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][19]))
                except:
                    text = "\n\033[1m\033[91mThe file\033[0m\n" + filename_history\
                           + "\n\033[1m\033[91mis corrupted. muLAn killed.\033[0m"
                    sys.exit(text)
            file.close()

        # Read on the chains
        if nb_chains > 0:
            for i in xrange(nb_chains):

                file = open(fnames_chains[i], 'r')
                for line in file:
                    params_model = line

                    if params_model[0] == '#':
                        continue

                    try:
                        samples_file['id'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][0]))
                        samples_file['t0'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][1]))
                        samples_file['u0'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][2]))
                        samples_file['tE'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][3]))
                        samples_file['rho'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][4]))
                        samples_file['gamma'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][5]))
                        samples_file['piEN'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][6]))
                        samples_file['piEE'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][7]))
                        samples_file['s'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][8]))
                        samples_file['q'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][9]))
                        samples_file['alpha'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][10]))
                        samples_file['dalpha'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][11]))
                        samples_file['ds'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][12]))
                        samples_file['chi2'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][13]))
                        samples_file['accrate'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][14]))
                        samples_file['date_save'].append(int(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][15]))
                        samples_file['time_save'].append(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][16])
                        samples_file['chi2/dof'].append(float(
                            [a for a in (params_model.split('\n')[0].split(' ')) if
                             (a != '')][17]))
                        samples_file['chain'].append(int(fnames_chains[i][-8:-4]))
                        samples_file['fullid'].append(-1)
                    except:
                        text = "\n\033[1m\033[91mThe file\033[0m\n" + "\033[1m\033[91m" + fnames_chains[i]\
                               + "\033[0m\n\033[1m\033[91mis corrupted. muLAn killed.\033[0m"
                        sys.exit(text)
                file.close()

        # Order the models
        chi2_min = np.min(samples_file['chi2'])
        samples_file.update(
            {'dchi2': samples_file['chi2'] - chi2_min})

        # Remove duplicates
        results = pd.DataFrame({})
        for key in samples_file:
            results[key] = samples_file[key]

        results_sorted = results.sort_values(['dchi2', 'fullid'], ascending=[1, 0]).drop_duplicates(
                subset=['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds', 'chi2'])

        # Give a unique ID to models
        id_start = np.max(results_sorted['fullid']) + 1
        if id_start == 0 : id_start = 1
        cond = results_sorted['fullid'] == -1
        results_sorted.loc[cond, 'fullid'] = id_start + np.arange(cond.sum())

        # Save new file in csv with exponential
        filename_history = filename_history[:-3] + 'csv'
        file = open(filename_history, 'w')
        format = '#{:},'.format('UniqueID')\
            + '{:},'.format('t0')\
            + '{:},'.format('u0')\
            + '{:},'.format('tE')\
            + '{:},'.format('rho')\
            + '{:},'.format('gamma')\
            + '{:},'.format('piEN')\
            + '{:},'.format('piEE')\
            + '{:},'.format('s')\
            + '{:},'.format('q')\
            + '{:},'.format('alpha')\
            + '{:},'.format('dalpha')\
            + '{:},'.format('ds')\
            + '{:},'.format('chi2')\
            + '{:},'.format('chi2/dof')\
            + '{:},'.format('accrate')\
            + '{:}'.format('chain')\
            + '\n'
        file.write(format)

        for i in xrange(len(results_sorted)):
            format = '{:},'.format(results_sorted['fullid'].values[i])\
                + '{:.10e},'.format(results_sorted['t0'].values[i])\
                + '{:.10e},'.format(results_sorted['u0'].values[i])\
                + '{:.10e},'.format(results_sorted['tE'].values[i])\
                + '{:.10e},'.format(results_sorted['rho'].values[i])\
                + '{:.10e},'.format(results_sorted['gamma'].values[i])\
                + '{:.10e},'.format(results_sorted['piEN'].values[i])\
                + '{:.10e},'.format(results_sorted['piEE'].values[i])\
                + '{:.10e},'.format(results_sorted['s'].values[i])\
                + '{:.10e},'.format(results_sorted['q'].values[i])\
                + '{:.10e},'.format(results_sorted['alpha'].values[i])\
                + '{:.10e},'.format(results_sorted['dalpha'].values[i])\
                + '{:.10e},'.format(results_sorted['ds'].values[i])\
                + '{:.10e},'.format(results_sorted['chi2'].values[i])\
                + '{:.10e},'.format(results_sorted['chi2/dof'].values[i])\
                + '{:.3f},'.format(results_sorted['accrate'].values[i])\
                + '{:}'.format(results_sorted['chain'].values[i])\
                + '\n'
            file.write(format)
        file.close()

# ----------------------------------------------------------------------
if (__name__ == "__main__"):
    pass

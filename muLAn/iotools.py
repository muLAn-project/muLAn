# -*-coding:Utf-8 -*
# ====================================================================
# Packages
# ====================================================================
import configparser as cp
import copy
import glob
import muLAn
import muLAn.packages.general_tools as gtools
import numpy as np
import os
import pandas as pd
import sys

class FitResults:
    """Class to read, save, and manipulate models tested during the fit.

    Args:
        parser (:obj:`configparser.ConfigParser`): options and configurations 
            for muLAn.
        run_id (str): Name of a muLAn archive (i.e., the name of a run).
            Default `None`.
        format (str): {`ascii` | `h5`}, default `ascii`. File format to load
            the MCMC results. 
    
    Attributes:
        samples (`pandas.DataFrame`): table of all the samples explored by the
            MCMC.

    """

    def __init__(self, parser, format='ascii', **kwargs):

        self.parser = parser

        # Load 
        if format=='ascii':
            self.load_aimc_from_file(parser, **kwargs)
        elif format=='h5':
            self.load()

    def load_aimc_from_file(self, cfgsetup, **kwargs):
        """Method to load model parameters from files created during MCMC.

        This method loads ASCII files created by the package EMCEE, after the
        end of an MCMC run. The sampler is assumed to be and AIMC.

        Args:
            parser (:obj:`configparser.ConfigParser`): options and configurations 
                for muLAn.
            run_id (str): Name of a muLAn archive (i.e., the name of a run).
                Default `None`.

        """

        # Identify filenames from MCMC
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')

        if 'run_id' in kwargs:
            fnames_chains = glob.glob(path + kwargs['run_id'] + "*-c*.txt")
            fnames_chains_exclude = glob.glob(path + kwargs['run_id'] + "*g*.txt")
        else:
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

            # Read on the chains
            if nb_chains > 0:
                for i in range(nb_chains):

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

            # Create a pandas.DataFrame to store the runs
            samples = pd.DataFrame(samples_file)
            samples['dchi2'] = samples['chi2'] - np.min(samples['chi2'])
            samples = samples.sort_values(['dchi2', 'fullid'], ascending=[1, 0])

            # Give a unique ID to models
            id_start = np.max(samples['fullid']) + 1
            if id_start == 0 : id_start = 1
            mask = samples['fullid'] == -1
            samples.loc[mask, 'fullid'] = id_start + np.arange(mask.sum())
            self.samples = samples

    def save(self, filename=None, format='h5', N=None):
        """Save MCMC samples in the specified format.

        Args:
            filename (str): file name for the output file.
            format (str, default 'h5'): {'ascii' | 'h5'}

        """
        if format == 'h5':
            if filename==None:
                fname = "{:s}-Fits.h5".format(self.parser.get('Controls', 'Archive'))
            else:
                fname = filename
            self.samples.to_hdf(fname, 'fits', mode='w')

        elif format == 'ascii':
            if filename==None:
                fname = "{:s}-Fits.csv".format(self.parser.get('Controls', 'Archive'))
            else:
                fname = filename

            if N == None:
                N = len(self.samples)

            # Save new file in csv with exponential
            file = open(fname, 'w')
            format = '#{:},'.format('UniqueID')\
                + '{:},'.format('dchi2')\
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

            for i in range(N):
                format = '{:},'.format(self.samples['fullid'].values[i])\
                    + '{:.3f},'.format(self.samples['dchi2'].values[i])\
                    + '{:.10e},'.format(self.samples['t0'].values[i])\
                    + '{:.10e},'.format(self.samples['u0'].values[i])\
                    + '{:.10e},'.format(self.samples['tE'].values[i])\
                    + '{:.10e},'.format(self.samples['rho'].values[i])\
                    + '{:.10e},'.format(self.samples['gamma'].values[i])\
                    + '{:.10e},'.format(self.samples['piEN'].values[i])\
                    + '{:.10e},'.format(self.samples['piEE'].values[i])\
                    + '{:.10e},'.format(self.samples['s'].values[i])\
                    + '{:.10e},'.format(self.samples['q'].values[i])\
                    + '{:.10e},'.format(self.samples['alpha'].values[i])\
                    + '{:.10e},'.format(self.samples['dalpha'].values[i])\
                    + '{:.10e},'.format(self.samples['ds'].values[i])\
                    + '{:.10e},'.format(self.samples['chi2'].values[i])\
                    + '{:.10e},'.format(self.samples['chi2/dof'].values[i])\
                    + '{:.3f},'.format(self.samples['accrate'].values[i])\
                    + '{:}'.format(self.samples['chain'].values[i])\
                    + '\n'
                file.write(format)
            file.close()


    def load(self, filename=None, format='h5'):
        """Save MCMC samples in the specified format.

        Args:
            fname (str): file name for the output file.
            format (str): currently, the only option is 'hdf5'.

        """
        if format == 'h5':
            if filename==None:
                fname = "{:s}-Fits.h5".format(self.parser.get('Controls', 'Archive'))
            else:
                fname = filename
            self.samples = pd.read_hdf(fname, 'fits')

    def remove_duplicates(self, inplace=False, **kwargs):
        """Create a table of MCMC samples without duplicates.

        Args:
            inplace (bool): default False. Replace self.samples if True, return
                the resulting table otherwise.

        """

        col = ['chi2', 't0', 'u0', 'tE', 'rho', 'gamma', 'piEE', 'piEN', 's',
                'q', 'alpha', 'dalpha', 'ds']
        if inplace:
            self.samples = self.samples.loc[
                    self.samples[col].round(12).drop_duplicates(subset=col).index]
        else:
            samples = self.samples.loc[
                    self.samples[col].round(12).drop_duplicates(subset=col).index]
            return samples


if (__name__ == "__main__"):
    pass


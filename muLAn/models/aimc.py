# -*-coding:Utf-8 -*

import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import copy
import emcee
import pickle
import glob
import shutil
import datetime
import importlib
import subprocess
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy import stats
from scipy import interpolate
from sklearn import linear_model
import muLAn.models as mulanmodels
import muLAn.packages.algebra as algebra


# Global variables to speed up the code
# -------------------------------------
global time_serie, fitp, constp, ndim
global fitpn, constpn, ndim, all_params
global instrument, algo, model_lib

# Ugly part of the code: some static arguments are loaded and declared as
# global to speed up the algorithms.
fname = 'args.h5'
data = pd.read_hdf(fname, 'data')
fitp = pd.read_hdf(fname, 'fit_params')
constp = pd.read_hdf(fname, 'const_params')

ndim = len(fitp)
fitpn = list(fitp.index)
constpn = list(constp.index)
instrument = np.unique(data['obs'])
algo = np.unique(data['model'])

all_params = fitp.to_dict()
all_params.update(constp.to_dict())

# Load library of models
model_lib = dict()
for i in range(algo.shape[0]):
    name = 'muLAn.models.{:s}'.format(algo[i])
    model_lib.update({algo[i]: importlib.import_module(name)})
# End of declaration of the global variables


def help():
    text = "AIMC - Affine Invariant MCMC."
    return text

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
            print(text)
        else:
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print(text)
# ----------------------------------------------------------------------
def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))


def logprior(param_model):
    p = 0
    if param_model['t0'] < 0:
        p = 1e12
    if param_model['rho'] < 0:
        p = 1e12
    if param_model['rho'] > 1.0:
        p = 1e12
    if param_model['tE'] < 1e-10:
        p = 1e12
    if param_model['q'] < 1e-9:
        p = 1e12
#    if param_model['q'] > 1.0:
#        p = 1e12
    if param_model['s'] < 1e-10:
        p = 1e12
    if param_model['s'] > 10:
        p = 1e12
    return p


def loglike(all_params, cfg):

    # For backwards compatibility (Deprecated)
    tb = all_params['tb']

    # Compute magnification 
    for j in range(len(instrument)):
        mask1 = data['obs'] == instrument[j]
        for i in range(algo.shape[0]):
            mask = (data['obs'] == instrument[j])\
                    & (data['model'] == algo[i])

            if mask.sum() > 0:
                epochs = data.loc[mask, 'dates'].values
                DsN = data.loc[mask, 'DsN'].values
                DsE = data.loc[mask, 'DsE'].values
                Ds = dict({'N': DsN, 'E': DsE})

                try:
                    kwargs_method = dict(parser.items(algo[i]))
                except:
                    kwargs_method = dict()

                mag = model_lib[algo[i]].magnifcalc(epochs, all_params, Ds=Ds, tb=tb, **kwargs_method)

                data.loc[mask,'amp'] = mag

        fs, fb = algebra.fsfbwsig(data[mask1], None, blending=True)
        data.loc[mask1,'fs'] = fs
        data.loc[mask1,'fb'] = fb

    data['flux_model'] = data['fs'] * data['amp'] + data['fb']
    chi2pp = np.power((data['flux']-data['flux_model'])/data['err_flux'], 2)
    chi2 = np.sum(chi2pp)
    result = - 0.5 * chi2
    
    return result

def logprob(theta, cfg):

    # Update the parameters
    for i in range(ndim):
        all_params[fitpn[i]] = theta[i]

    # Evaluate log likelihood
    return loglike(all_params, cfg) + logprior(all_params)


# ----------------------------------------------------------------------
def ini_chains_gene(cfg):

    nwalkers = cfg.getint('AIMC', 'walkers')

    result = []
    j = 0
    while(j<nwalkers):
        table = np.array([])
        for i in range(ndim):
            l = cfg.get('Modelling', fitpn[i])
            a = abs(float(l.split(',')[1]))
            b = abs(float(l.split(',')[2]))
            c = float(l.split(',')[3])
            aa = c-a
            bb = c+b
            x = (bb - aa) * np.random.random_sample() + aa
            table = np.append(table, x)
        result.append(table)
        j+=1

    return result


def search(**kwargs):

    # Declare global variables
    global ndim, fitp

    # Parse the arguments
    if 'cfgsetup' in kwargs: cfgsetup = kwargs['cfgsetup']

    # Create alias for paths
    path_event = cfgsetup.get('FullPaths', 'Event')
    path_mcmc = f"{path_event}/{cfgsetup.get('RelativePaths', 'Chains')}"
    archive_name = cfgsetup.get('Controls', 'Archive')

    # Extract flags from parser
    flag_resume = cfgsetup.getboolean('AIMC', 'resume')

    # Extract MCMC settings
    nwalkers = cfgsetup.getint('AIMC', 'walkers')
    length = cfgsetup.getint('AIMC', 'length')
    ncpu = cfgsetup.getint('AIMC', 'cpu')

    # Create a file to check if user wants to stop MCMC
    fn_lock = cfgsetup.get('FullPaths', 'Event') + '.lock'
    if not os.path.exists(fn_lock): open(fn_lock, 'w').close()

    if not flag_resume:

        # Initialize folders tree
        shutil.rmtree(path_mcmc)
        if not os.path.exists(path_mcmc): os.makedirs(path_mcmc)

        # Set up the backend
        backend = emcee.backends.HDFBackend(f'{path_mcmc}/mcmc_sampler.h5')
        backend.reset(nwalkers, ndim)

        # Initialize chains
        pos = ini_chains_gene(cfgsetup)

        with Pool(processes=ncpu) as pool:
            
            # Initialize the sampler
            sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob, args=[cfgsetup], backend=backend, pool=pool)

            # Now we start sampling
            for sample in sampler.sample(pos, iterations=length, progress=True):

                # Check if user wants to stop the MCMC
                if not os.path.exists(fn_lock): break
    else:
        backend = emcee.backends.HDFBackend(f'{path_mcmc}/mcmc_sampler.h5')
        print("Initial size: {0}".format(backend.iteration))
        with Pool(processes=ncpu) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob, args=[cfgsetup], backend=backend, pool=pool)
            sampler.run_mcmc(None, length, progress=True)

    if os.path.exists(fn_lock): os.remove(fn_lock)
    else: sys.exit("\nProcess stopped by the user.\n")


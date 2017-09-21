# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# Plot correlations between the parameters
# ----------------------------------------------------------------------
#   External libraries
# ----------------------------------------------------------------------
import sys
import os

# Full path of this file
full_path_here = os.path.realpath(__file__)
text = full_path_here.split('/')
a = ''
i = 0
while i < len(text) - 1:
    a = a + text[i] + '/'
    i = i + 1
full_path = a

#filename = full_path + '../' + '.pythonexternallibpath'
#file = open(filename, 'r')
#for line in file:
#    path_lib_ext = line
#file.close()
#if path_lib_ext != 'None':
#    sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#   Standard packages
# ----------------------------------------------------------------------
import os
import glob
import getdist as getdist
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
import warnings
import sys
import copy
import cmath
import shutil
import math
import emcee
# import pylab
import pickle
import pylab
import zipfile
import datetime
import subprocess
import numpy as np
import pandas as pd
# import bokeh as bok
from scipy import stats
import ConfigParser as cp
from astropy.time import Time
from PyAstronomy import pyasl
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.ticker import FixedLocator, FormatStrFormatter
# ----------------------------------------------------------------------
#   Non-standard packages
# ----------------------------------------------------------------------
# import models.ephemeris as ephemeris
# import models.esblparall as esblparall
# import packages.plotconfig as plotconfig
# import models.esblparallax as esblparallax
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
# -----------------------------------------------------------------------------
class cufloat:
    """Uncertainty"""
    def __init__(self,x,sp,sm):
        self.x = x
        self.sp = sp
        self.sm = sm

    def output(self):
        xp = ufloat(self.x, self.sp)
        xm = ufloat(self.x, self.sm)

        xp_out = '{:}'.format(xp)
        xm_out = '{:}'.format(xm)

        xp_split = xp_out.split("+/-")
        xm_split = xm_out.split("+/-")

        if(xp_split[0][0] == "("):
            flag_e = 1
            tempp = xp_out.split("e")
            tempm = xm_out.split("e")
            exposant = tempp[1]

            xp_out = tempp[0][1:-1]
            xm_out = tempm[0][1:-1]

            xp_split = xp_out.split("+/-")
            xm_split = xm_out.split("+/-")
        else: flag_e = 0


        nb_p = len(xp_split[0])
        nb_m = len(xm_split[0])
        if(nb_p>=nb_m): x_cen = xp_split[0]
        else: x_cen = xm_split[0]

        if(flag_e==0):
            final = x_cen + "+/-[" + xp_split[1] + "/" + xm_split[1] + "]"
        if(flag_e==1):
            final = "(" + x_cen + "+/-[" + xp_split[1] + "/" + xm_split[1] + "])x10^[" + exposant + "]"
#            final = "(" + x_cen + "+/-[" + xp_split[1] + "/" + xm_split[1] + "])e[" + exposant + "]"

        return final
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
def help():
    text = "Plot the correlations between the parameters using the external package GetDist."
    return text


# ----------------------------------------------------------------------
def plot(cfgsetup=False, models=False, model_param=False, time_serie=False,
         obs_properties=False, options=False, interpol_method=False):

    # -----------------------------------------------------------------------------
    #   Annonce
    # -----------------------------------------------------------------------------
    # text = "\n"
    # communicate(cfgsetup, 1, text, opts=False)
    # text = "Plotting the correlations"
    # options = [printoption.reverse]
    # communicate(cfgsetup, 1, text, opts=options)

    # -----------------------------------------------------------------------------
    #   Load parameters
    # -----------------------------------------------------------------------------
    # Parameter to be fitted / Nuisance parameters
    params = {
        't0' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    't0').split(',')]),\
        'u0' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'u0').split(',')]),\
        'tE' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'tE').split(',')]),\
        'rho' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'rho').split(',')]),\
        'gamma' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'gamma').split(',')]),\
        'piEE' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'piEE').split(',')]),\
        'piEN' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'piEN').split(',')]),\
        's' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    's').split(',')]),\
        'q' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'q').split(',')]),\
        'alpha' : np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'alpha').split(',')]),\
        'dalpha': np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'dalpha').split(',')]),\
        'ds': np.array([a.strip() for a in cfgsetup.get('Modelling',
                    'ds').split(',')])\
        }

    fitted_param = dict()
    result = np.array([])
    if params['t0'][0] == "fit":
        fitted_param.update({'t0': params['t0'][3].astype(np.float64)})
        result = np.append(result, fitted_param['t0'])
    if params['u0'][0] == "fit":
        fitted_param.update({'u0': params['u0'][3].astype(np.float64)})
        result = np.append(result, fitted_param['u0'])
    if params['tE'][0] == "fit":
        fitted_param.update({'tE': params['tE'][3].astype(np.float64)})
        result = np.append(result, fitted_param['tE'])
    if params['rho'][0] == "fit":
        fitted_param.update({'rho': params['rho'][3].astype(np.float64)})
        result = np.append(result, fitted_param['rho'])
    if params['gamma'][0] == "fit":
        fitted_param.update({'gamma': params['gamma'][3].astype(np.float64)})
        result = np.append(result, fitted_param['gamma'])
    if params['piEE'][0] == "fit":
        fitted_param.update({'piEE': params['piEE'][3].astype(np.float64)})
        result = np.append(result, fitted_param['piEE'])
    if params['piEN'][0] == "fit":
        fitted_param.update({'piEN': params['piEN'][3].astype(np.float64)})
        result = np.append(result, fitted_param['piEN'])
    if params['s'][0] == "fit":
        fitted_param.update({'s': params['s'][3].astype(np.float64)})
        result = np.append(result, fitted_param['s'])
    if params['q'][0] == "fit":
        fitted_param.update({'q': params['q'][3].astype(np.float64)})
        result = np.append(result, fitted_param['q'])
    if params['alpha'][0] == "fit":
        fitted_param.update({'alpha': params['alpha'][3].astype(np.float64)})
        result = np.append(result, fitted_param['alpha'])
    if params['dalpha'][0] == "fit":
        fitted_param.update({'dalpha': params['dalpha'][3].astype(np.float64)})
        result = np.append(result, fitted_param['dalpha'])
    if params['ds'][0] == "fit":
        fitted_param.update({'ds': params['ds'][3].astype(np.float64)})
        result = np.append(result, fitted_param['ds'])

    # -----------------------------------------------------------------------------
    #   Load the chain files
    # -----------------------------------------------------------------------------

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
    if nb_chains > 0:
        for i in xrange(nb_chains):
            file = open(fnames_chains[i], 'r')
            samples_file = dict({'chi2': [], 't0': [], 'u0': [], 'tE': [], 'rho': [], 'gamma': [], 'piEE': [], 'piEN': [], 's': [], 'q': [], 'alpha': [], 'dalpha': [], 'ds': [], 'chain': [], 'fullid': [], 'date_save': [], 'time_save': [], 'id': [], 'accrate': [], 'chi2/dof': []})

            for line in file:
                params_model = line

                if params_model[0] == '#':
                    continue

                samples_file['id'].append(int([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][0]))
                samples_file['t0'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][1]))
                samples_file['u0'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][2]))
                samples_file['tE'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][3]))
                samples_file['rho'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][4]))
                samples_file['gamma'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][5]))
                samples_file['piEN'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][6]))
                samples_file['piEE'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][7]))
                samples_file['s'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][8]))
                samples_file['q'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][9]))
                samples_file['alpha'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][10]))
                samples_file['dalpha'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][11]))
                samples_file['ds'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][12]))
                samples_file['chi2'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][13]))
                samples_file['accrate'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][14]))
                samples_file['date_save'].append(int([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][15]))
                samples_file['time_save'].append([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][16])
                samples_file['chi2/dof'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][17]))
                samples_file['chain'].append(int(fnames_chains[i][-8:-4]))
                samples_file['fullid'].append(-1)

            file.close()

            if i==0:
                chains = pd.DataFrame({})
                for key in samples_file:
                    chains[key] = samples_file[key]
            else:
                chains_tmp = pd.DataFrame({})
                for key in samples_file:
                    chains_tmp[key] = samples_file[key]
                chains = chains.append(chains_tmp, ignore_index=True)
                del chains_tmp

    chains['reject'] = np.full(len(chains), 0, dtype='i4')

    # -----------------------------------------------------------------------------
    #   Write files for GetDist package
    # -----------------------------------------------------------------------------
    path_formatted = path + "formatted-chains/"
    if os.path.exists(path_formatted):
        shutil.rmtree(path_formatted)

    if not os.path.exists(path_formatted):
        os.makedirs(path_formatted)

    order = np.array(['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds'])
    names_leg = np.array(['\mathsf{t_0}', '\mathsf{u_0}', '\mathsf{t_E}', '\\rho', '\Gamma', '\pi_\mathsf{EN}', '\pi_\mathsf{EE}', '\mathsf{s}', '\mathsf{q}', '\\alpha', '\mathrm{d}\\alpha/\mathrm{d}t', '\mathrm{d}s/\mathrm{d}t'])

    nb_params = len(fitted_param)
    names = np.array([key for key in fitted_param])
    names2 = np.array([])
    names2_leg = np.array([])
    for i in xrange(len(order)):
        cond = np.where(order[i] == names)[0]
        if len(cond) > 0:
            names2 = np.append(names2, order[i])
            names2_leg = np.append(names2_leg, names_leg[i])

    names = copy.deepcopy(names2)
    names_leg = copy.deepcopy(names2_leg)
    del names2
    del names2_leg

    text = ""
    for i in xrange(len(names)):
        text = text + "{} {:.3f} {:.3f}\n".format(names[i], -1e10, 1e10)

    filename = "{:}.ranges".format(path_formatted + cfgsetup.get('Controls', 'Archive'))
    file = open(filename, "w")
    file.write(text)
    file.close()

    text = ""
    for i in xrange(len(names)):
        text = text + "{:} {:}\n".format(names[i], names_leg[i])

    filename = "{:}.paramnames".format(path_formatted + cfgsetup.get('Controls', 'Archive'))
    file = open(filename, "w")
    file.write(text)
    file.close()

    if nb_chains > 0:
        for i in xrange(nb_chains):

            cond = chains.chain == i
            chain_experiment = chains[cond]

            text = ""
            for j in xrange(len(chain_experiment)):
                if chain_experiment['chi2'].values[j] != np.inf:
                    text = text + "1 {} ".format(repr(0.5*chain_experiment['chi2'].values[j]))
                    for k in xrange(nb_params):
                        text = text + "{:} ".format(repr(chain_experiment[names[k]].values[j]))
                    text = text + "\n"

            filename = "{:}_{:d}.txt".format(path_formatted + cfgsetup.get('Controls', 'Archive'), i+1)
            file = open(filename, "w")
            file.write(text)
            file.close()

    # -----------------------------------------------------------------------------
    #   Load samples files in the GetDist conventions
    # -----------------------------------------------------------------------------
    # path_formatted = path + "formatted-chains-clean-renamed/"
    from getdist import mcsamples, plots
    # Load
    print "      The output is now redirected to:\n      {:}.".format(path_formatted + "output_getdist.log")
    filename = "{:}{:}".format(path_formatted, cfgsetup.get('Controls', 'Archive'))
    terminal_output = sys.stdout
    terminalerr_output = sys.stderr
    sys.stdout = open(path_formatted + "corr02_output.txt", "w")
    sys.stderr = open(path_formatted + "corr02_errors.log", "w")
    print "Sample files loaded :"
    samples = mcsamples.loadMCSamples(filename, settings={'ignore_rows': 0.0})
    filename_conv = path + "formatted-chains/corr02_getdist_output.txt"
    samples.getConvergeTests(test_confidence=0.95, writeDataToFile=True, what=['MeanVar', 'GelmanRubin', 'SplitTest', 'RafteryLewis', 'CorrLengths'], filename=filename_conv, feedback=True)
    # Plot object
    g = plots.getSubplotPlotter()
    # Scatter plot
    g.triangle_plot(samples, filled=False)
    # Best model
    print samples.getLikeStats()
    # print samples.getLikeStats().names[0].bestfit_sample
    # Add the samples on the plots and the best-fitting parameters
    chains['dchi2'] = chains['chi2'] - np.min(chains['chi2'].values)
    # dchi2_list = np.array([999999999, 100, 50, 10, 0])
    dchi2_list = np.array([1e12, 18.4, 11.8, 6.17, 2.30, 0.0])
    # colors = ['lightsteelblue', 'deepskyblue', 'limegreen', 'red']
    spectral_map = plt.get_cmap('Spectral')
    colors = ['black', spectral_map(255.0/255.0), spectral_map(170.0/255.0), spectral_map(85.0/255.0), spectral_map(0.0)]
    # colors = np.array(colors)/255.0
    alpha = [0.5, 0.9, 0.9, 0.9, 0.9]
    s_list = [1, 5, 5, 5, 5]
    for id_dchi2 in xrange(len(dchi2_list)-1):
        cond = chains.reject == 1
        chains.loc[cond,'reject'] = 0
        chains.loc[((chains.dchi2 < dchi2_list[id_dchi2+1]) | (chains.dchi2 >= dchi2_list[id_dchi2])), 'reject'] = 1
        cond = chains.reject == 0
        for j in xrange(len(g.subplots) - 1):
            for i in xrange(len(g.subplots)-j-1):
                g.subplots[i+j+1, j].scatter(chains[cond][names[j]].values, chains[cond][np.roll(names, -j-1)[i]].values, s=s_list[id_dchi2], facecolors=colors[id_dchi2], marker='o', alpha=alpha[id_dchi2], linewidths=0, zorder=-100+id_dchi2)
                g.subplots[i+j+1, j].axvline(samples.getLikeStats().names[j].bestfit_sample, color='black', dashes=(4, 3), lw=0.5, zorder=1)
                g.subplots[i+j+1, j].axhline(samples.getLikeStats().names[i+j+1].bestfit_sample, color='black', dashes=(4, 3), lw=0.5, zorder=1)

    path_plots = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Plots')
    filename = "{:s}{:s}-corr.png".format(path_plots, cfgsetup.get('Controls', 'Archive'))
    g.export(filename)

    sys.stdout.close()
    sys.stdout = terminal_output

    sys.stderr.close()
    sys.stderr = terminalerr_output
    print "      The output is now redirected to the terminal."

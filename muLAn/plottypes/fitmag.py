# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# Routine to plot the result of the MCMC, in magnitude.
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

# filename = full_path + '../' + '.pythonexternallibpath'
# file = open(filename, 'r')
# for line in file:
#     path_lib_ext = line
# file.close()
# if path_lib_ext != 'None':
#     sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#  Standard packages
# ----------------------------------------------------------------------
import os
import glob
import sys
import copy
import cmath
# import math
import emcee
# import pylab
import pickle
import pylab
import zipfile
import datetime
from scipy import interpolate
import subprocess
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.optimize import fsolve
import pandas as pd
import bokeh.layouts as blyt
import bokeh.plotting as bplt
from bokeh.models import HoverTool, TapTool, ColumnDataSource, OpenURL
from bokeh.models.widgets import DateFormatter, NumberFormatter, DataTable, \
    TableColumn
import bokeh.io as io
from scipy import stats
import ConfigParser as cp
from astropy.time import Time
from PyAstronomy import pyasl
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.ticker import FixedLocator, FormatStrFormatter
# ----------------------------------------------------------------------
#  Non-standard packages
# ----------------------------------------------------------------------
import muLAn.models.ephemeris as ephemeris


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
    text = "Plot the light curve of a previously modelled event."
    return text
# ----------------------------------------------------------------------

def bash_command(text):
    proc = subprocess.Popen(text, shell=True, executable="/bin/bash")
    proc.wait()


# ----------------------------------------------------------------------
def unpack_options(cfgsetup, level0, level1, sep=','):
    options = [a.strip() for a in cfgsetup.get(level0, level1).split(sep)]
    del a, cfgsetup, level0, level1
    return options


# ----------------------------------------------------------------------
def fsfb(time_serie, cond, blending=True):

    #blending = True

    x = np.atleast_2d(time_serie['amp'][cond]).T
    y = np.atleast_2d(time_serie['flux'][cond]).T

    regr = linear_model.LinearRegression(fit_intercept=blending)
    regr.fit(x, y)
    fs = regr.coef_[0][0]
    # fb = regr.intercept_[0]
    if blending:
        fb = regr.intercept_[0]
    else:
        fb = 0.0

    return fs, fb

# ----------------------------------------------------------------------
def critic_roots(s, q, phi):
    """Sample of the critic curve. The convention is :
    - the heaviest body (mass m1) is the origin;
    - the lightest body (mass m2) is at (-s, 0).

    Arguments:
    s -- the binary separation;
    q -- the lens mass ratio q = m2/m1;
    phi -- the sample parameter in [0;2*pi].

    Returns:
    result -- numpy array of the complex roots.
    """

    coefs = [1, 2 * s, s ** 2 - np.exp(1j * phi),
             -2 * s * np.exp(1j * phi) / (1 + q),
             -(s ** 2 * np.exp(1j * phi) / (1 + q))]
    result = np.roots(coefs)
    del coefs
    return result


# ----------------------------------------------------------------------
# Levi-Civita coefficient
def epsilon(i, j, k):
    if (i == j or i == k or j == k):
        e = 0
    else:
        if (i == 1 and j == 2 and k == 3):
            e = 1
        if (i == 3 and j == 1 and k == 2):
            e = 1
        if (i == 2 and j == 3 and k == 1):
            e = 1
        if (i == 1 and j == 3 and k == 2):
            e = -1
        if (i == 3 and j == 2 and k == 1):
            e = -1
        if (i == 2 and j == 1 and k == 3):
            e = -1
    return e


#
# Projection onto the sky
def onSky(m_hat, n_hat, u):
    x = np.array(
        [epsilon(i + 1, j + 1, k + 1) * u[i] * n_hat[j] * m_hat[k] for i
         in xrange(3) for j in xrange(3) for k in xrange(3)]).sum()
    u_proj = (1.0 / np.sqrt(1 - ((n_hat * m_hat).sum()) ** 2)) \
             * np.array([x,
                         (n_hat * u).sum() - (n_hat * m_hat).sum() * (
                             u * m_hat).sum(), 0])
    return u_proj


#
# Projection onto the sky
def normalize(u):
    return u / np.sqrt((u * u).sum())


def angle_between(v1, v2):
    v1_u = normalize(v1)
    v2_u = normalize(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


#
#
def vectoriel(u, v):
    x = u[1] * v[2] - u[2] * v[1]
    y = u[2] * v[0] - u[0] * v[2]
    z = u[0] * v[1] - u[1] * v[0]
    return np.array([x, y, z])


def boussole(EarthSunFile=False, EarthSatelliteFile=False, cfg=False, \
             t_D_xy=False):
    # Value of the origin of the developments
    t0par = cfg.getfloat('Modelling', 'tp')

    # Coordinates conversion of the event from the Equatorial frame to the Ecliptic frame
    c_icrs = SkyCoord(ra=cfg.get('EventDescription', 'RA'),
                      dec=cfg.get('EventDescription', 'DEC'),
                      frame='icrs')
    l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
    b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

    # Vector Earth --> Sun in Ecliptic frame (cartesian coordinates).
    # ------------------------------------------------------------------
    format = {'names': ('dates', 'x', 'y', 'z', 'vx', 'vy', 'vz'), \
              'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}

    temp = np.loadtxt(EarthSunFile, usecols=(0, 5, 6, 7, 8, 9, 10),
                      dtype=format, unpack=False)
    EarthSun = pd.DataFrame(temp)

    del temp

    # Time conversion: TDB->TCG->HJD
    temp = EarthSun['dates'] - 2400000.0
    flag_clem = 0
    if flag_clem:
        EarthSun['hjd'] = np.array(
            [pyasl.helio_jd(tc, c_icrs.ra.degree, c_icrs.dec.degree) for
             tc in
             Time(temp, format='mjd', scale='tdb').tcg.value]) - 50000.0
    else:
        EarthSun['hjd'] = np.array([pyasl.helio_jd(tc, l, b) for tc in
                                    Time(temp, format='mjd',
                                         scale='tdb').tcg.value]) - 50000.0

    del temp

    # Vector Earth --> Satellite in Ecliptic frame (cartesian coordinates).
    # ------------------------------------------------------------------
    format = {'names': ('dates', 'x', 'y', 'z'), \
              'formats': ('f8', 'f8', 'f8', 'f8')}

    temp = np.loadtxt(EarthSatelliteFile, usecols=(0, 5, 6, 7),
                      dtype=format, unpack=False)
    EarthSat = pd.DataFrame(temp)
    del temp

    # Time conversion: TDB->TCG->HJD
    temp = EarthSun['dates'] - 2400000.0
    flag_clem = 0
    if flag_clem:
        EarthSat['hjd'] = np.array(
            [pyasl.helio_jd(tc, c_icrs.ra.degree, c_icrs.dec.degree) for
             tc in
             Time(temp, format='mjd', scale='tdb').tcg.value]) - 50000.0
    else:
        EarthSat['hjd'] = np.array([pyasl.helio_jd(tc, l, b) for tc in
                                    Time(temp, format='mjd',
                                         scale='tdb').tcg.value]) - 50000.0

    del temp

    # Vector Earth --> Sun and velocity( Earth --> Sun ) at t0par
    sp = np.array(
        [interp1d(EarthSun['hjd'], EarthSun['x'], kind='linear')(t0par), \
         interp1d(EarthSun['hjd'], EarthSun['y'], kind='linear')(t0par), \
         interp1d(EarthSun['hjd'], EarthSun['z'], kind='linear')(
             t0par)])

    vp = np.array([interp1d(EarthSun['hjd'], EarthSun['vx'],
                            kind='linear')(t0par),
                   interp1d(EarthSun['hjd'], EarthSun['vy'],
                            kind='linear')(t0par),
                   interp1d(EarthSun['hjd'], EarthSun['vz'],
                            kind='linear')(t0par)])

    # Ecliptic frame [gamma, y, nord], cartesian coordinates
    n_hat = np.array([0, 0, 1])
    m_hat = np.array([np.cos(np.radians(b)) * np.cos(np.radians(l)), \
                      np.cos(np.radians(b)) * np.sin(np.radians(l)), \
                      np.sin(np.radians(b))])

    # Sky ref. frame [East, North projected, microlens]
    # Cartesian coordinates
    # Parallax correction from Earth
    delta_pos = np.array([])
    delta_pos_proj = np.array([])
    pos_proj = np.array([])
    for t in xrange(len(EarthSun)):
        pos = np.array(
            [EarthSun['x'][t], EarthSun['y'][t], EarthSun['z'][t]])
        delta_pos_temp = pos - (EarthSun['hjd'][t] - t0par) * vp - sp
        delta_pos_proj = np.append(delta_pos_proj,
                                   onSky(m_hat, n_hat, delta_pos_temp))
        pos_proj = np.append(pos_proj, onSky(m_hat, n_hat, pos))
        delta_pos = np.append(delta_pos, delta_pos_temp)

    delta_pos = np.reshape(delta_pos, (delta_pos.shape[0] / 3, 3))
    pos_proj = np.reshape(pos_proj, (pos_proj.shape[0] / 3, 3))
    delta_pos_proj = np.reshape(delta_pos_proj,
                                (delta_pos_proj.shape[0] / 3, 3))

    EarthSun['xproj'] = pos_proj.T[0]
    EarthSun['yproj'] = pos_proj.T[1]
    EarthSun['zproj'] = pos_proj.T[2]

    EarthSun['deltaxproj'] = delta_pos_proj.T[0]
    EarthSun['deltayproj'] = delta_pos_proj.T[1]
    EarthSun['deltazproj'] = delta_pos_proj.T[2]

    EarthSun['deltax'] = delta_pos.T[0]
    EarthSun['deltay'] = delta_pos.T[1]
    EarthSun['deltaz'] = delta_pos.T[2]

    # Correction due to the Satellite + parallax
    delta_pos = np.array([])
    delta_pos_proj = np.array([])
    pos_proj = np.array([])
    for t in xrange(len(EarthSat)):
        pos = np.array([EarthSun['x'][t] - EarthSat['x'][t],
                        EarthSun['y'][t] - EarthSat['y'][t],
                        EarthSun['z'][t] - EarthSat['z'][t]])
        delta_pos_temp = pos - (EarthSat['hjd'][t] - t0par) * vp - sp
        delta_pos_proj = np.append(delta_pos_proj,
                                   onSky(m_hat, n_hat, delta_pos_temp))
        pos_proj = np.append(pos_proj, onSky(m_hat, n_hat, pos))
        delta_pos = np.append(delta_pos, delta_pos_temp)

    delta_pos = np.reshape(delta_pos, (delta_pos.shape[0] / 3, 3))
    pos_proj = np.reshape(pos_proj, (pos_proj.shape[0] / 3, 3))
    delta_pos_proj = np.reshape(delta_pos_proj,
                                (delta_pos_proj.shape[0] / 3, 3))

    EarthSat['xproj'] = pos_proj.T[0]
    EarthSat['yproj'] = pos_proj.T[1]
    EarthSat['zproj'] = pos_proj.T[2]

    EarthSat['deltaxproj'] = delta_pos_proj.T[0]
    EarthSat['deltayproj'] = delta_pos_proj.T[1]
    EarthSat['deltazproj'] = delta_pos_proj.T[2]

    EarthSat['deltax'] = delta_pos.T[0]
    EarthSat['deltay'] = delta_pos.T[1]
    EarthSat['deltaz'] = delta_pos.T[2]

    if t_D_xy != False:
        D_ecl = np.array([interp1d(EarthSat['hjd'], EarthSat['x'],
                                   kind='linear')(t_D_xy), \
                          interp1d(EarthSat['hjd'], EarthSat['y'],
                                   kind='linear')(t_D_xy), \
                          interp1d(EarthSat['hjd'], EarthSat['z'],
                                   kind='linear')(t_D_xy)])
        D_enm = onSky(m_hat, n_hat, D_ecl)


    # print D_enm
    return EarthSun, EarthSat, D_enm


# ----------------------------------------------------------------------
#    Functions used to visualise DMCMC results
# ----------------------------------------------------------------------
def plot(cfgsetup=False, models=False, model_param=False, time_serie=False, \
         obs_properties=False, options=False, interpol_method=False):

    #    Initialisation of parameters
    # ------------------------------------------------------------------
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
        'dalpha': np.array([a.strip() for a in cfgsetup.get('Modelling', 'alpha').split(',')]),\
        'ds': np.array([a.strip() for a in cfgsetup.get('Modelling', 'alpha').split(',')])\
        }

    flag_fix_gamma = 1
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
        flag_fix_gamma = 0
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

    nb_param_fit = len(fitted_param)

    #    Initialisation
    # ------------------------------------------------------------------
    path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths',
                                                             'Chains')

    fnames_chains = glob.glob(
        path + cfgsetup.get('Controls', 'Archive') + "*-c*.txt")
    fnames_chains_exclude = glob.glob(
        path + cfgsetup.get('Controls', 'Archive') + "*g*.txt")

    temp = []
    for a in fnames_chains:
        if (a in fnames_chains_exclude) == False:
            temp.append(a)

    fnames_chains = copy.deepcopy(temp)
    del temp, fnames_chains_exclude

    nb_chains = len(fnames_chains)

    samples_file = dict(
        {'chi2': [], 't0': [], 'u0': [], 'tE': [], 'rho': [], \
         'gamma': [], 'piEE': [], 'piEN': [], 's': [], 'q': [], \
         'alpha': [], 'dalpha': [], 'ds': [], 'chain': [], 'fullid': [], 'chi2': [], 'chi2/dof': [],\
         'date_save': [], 'time_save': [], 'id': [], 'accrate': []})

    # filename_history = cfgsetup.get('FullPaths', 'Event') \
    #                    + cfgsetup.get('RelativePaths', 'ModelsHistory') \
    #                    + 'ModelsHistory.txt'

    filename_history = cfgsetup.get('FullPaths', 'Event') \
                       + cfgsetup.get('RelativePaths', 'ModelsHistory') \
                       + cfgsetup.get('Controls', 'Archive') \
                       + '-ModelsSummary.csv'

    flag_fix = 0
    labels = ['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds']
    for lab in labels:
        if unpack_options(cfgsetup, 'Modelling', lab)[0]!='fix':
            flag_fix = 1

    if os.path.exists(filename_history) & flag_fix:
        file = open(filename_history, 'r')
        for line in file:
            params_model = line

            if params_model[0] == '#':
                continue

            samples_file['fullid'].append(int(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][0]))
            samples_file['t0'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][1]))
            samples_file['u0'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][2]))
            samples_file['tE'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][3]))
            samples_file['rho'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][4]))
            samples_file['gamma'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][5]))
            samples_file['piEN'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][6]))
            samples_file['piEE'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][7]))
            samples_file['s'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][8]))
            samples_file['q'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][9]))
            samples_file['alpha'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][10]))
            samples_file['dalpha'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][11]))
            samples_file['ds'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][12]))
            samples_file['chi2'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][13]))
            samples_file['chi2/dof'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][14]))
            samples_file['accrate'].append(float(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][15]))
            samples_file['chain'].append(int(
                [a for a in (params_model.split('\n')[0].split(',')) if
                 (a != '')][16]))
        file.close()
    elif flag_fix:
        # Read on the chains
        if nb_chains > 0:
            for i in xrange(nb_chains):

                file = open(fnames_chains[i], 'r')
                for line in file:
                    params_model = line

                    if params_model[0] == '#':
                        continue

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
                    samples_file['time_save'].append(int(
                        [a for a in (params_model.split('\n')[0].split(' ')) if
                         (a != '')][16]))
                    samples_file['chi2/dof'].append(float(
                        [a for a in (params_model.split('\n')[0].split(',')) if
                         (a != '')][17]))
                    samples_file['chain'].append(int(fnames_chains[i][-8:-4]))
                    samples_file['fullid'].append(-1)
                file.close()

    # TO BE REMOVE: PB with negative rho.
    for ii in xrange(len(samples_file['rho'])):
        if samples_file['rho'][ii] < 0:
            samples_file['rho'][ii] = 0.000001
    # ------------------------------------

    # Best model
    # ------------------------------------------------------------------
    rang_2plot = [0]
    if flag_fix:
        models2plot = unpack_options(cfgsetup, 'Plotting', 'Models')
        if len(models2plot)==1:
            models2plot = unpack_options(cfgsetup, 'Plotting', 'Models', sep='-')
            if len(models2plot) == 2:
                first = int(models2plot[0])
                second = int(models2plot[1])
                nb_local = second - first + 1
                models2plot = np.linspace(first, second, nb_local, endpoint=True, dtype='i4')
                models2plot = ['{:d}'.format(a) for a in models2plot]
        rang_2plot = []
        if (models2plot != ['']) & (os.path.exists(filename_history)):
            chi2_min = np.min(samples_file['chi2'])
            samples_file.update(
                {'dchi2': samples_file['chi2'] - chi2_min})
            rang_best_model = np.where(samples_file['dchi2'] == 0)[0][0]
            for mod in models2plot:
                if mod != '':
                    try:
                        rang_2plot.append(
                            np.where(np.array(samples_file['fullid']) == int(mod))[
                                0][0])
                    except:
                        sys.exit('Cannot find the models you ask me to plot.')
        else:
            chi2_min = np.min(samples_file['chi2'])
            samples_file.update(
                {'dchi2': samples_file['chi2'] - chi2_min})
            rang_best_model = np.where(samples_file['dchi2'] == 0)[0][0]
            rang_2plot = [rang_best_model]
    else:
        rang_best_model = 0

    # Plots
    for idmod in xrange(len(rang_2plot)):
        if flag_fix:
            best_model = dict({})
            best_model.update({'t0': samples_file['t0'][rang_2plot[idmod]]})
            best_model.update({'u0': samples_file['u0'][rang_2plot[idmod]]})
            best_model.update({'tE': samples_file['tE'][rang_2plot[idmod]]})
            best_model.update({'rho': samples_file['rho'][rang_2plot[idmod]]})
            best_model.update({'gamma': samples_file['gamma'][rang_2plot[idmod]]})
            best_model.update({'piEE': samples_file['piEE'][rang_2plot[idmod]]})
            best_model.update({'piEN': samples_file['piEN'][rang_2plot[idmod]]})
            # best_model.update({'piE': np.sqrt(np.power(samples_file['piEN'][rang_2plot[idmod]],2) + np.power(samples_file['piEN'][rang_2plot[idmod]],2))})
            best_model.update({'s': samples_file['s'][rang_2plot[idmod]]})
            best_model.update({'q': samples_file['q'][rang_2plot[idmod]]})
            best_model.update({'alpha': samples_file['alpha'][rang_2plot[idmod]]})
            best_model.update({'dalpha': samples_file['dalpha'][rang_2plot[idmod]]})
            best_model.update({'ds': samples_file['ds'][rang_2plot[idmod]]})
            best_model.update(
                {'chi2': samples_file['chi2'][rang_2plot[idmod]]})
            best_model.update(
                {'chi2/dof': samples_file['chi2/dof'][rang_2plot[idmod]]})
            # best_model.update({'id': samples_file['id'][rang_2plot[idmod]]})
            best_model.update({'chain': samples_file['chain'][rang_2plot[idmod]]})
            # best_model.update(
            #     {'date_save': samples_file['date_save'][rang_2plot[idmod]]})
            # best_model.update(
            #     {'time_save': samples_file['time_save'][rang_2plot[idmod]]})
            best_model.update(
                {'accrate': samples_file['accrate'][rang_2plot[idmod]]})
            best_model.update(
                {'fullid': samples_file['fullid'][rang_2plot[idmod]]})
        else:
            best_model = dict({})
            #samples_file = dict()
            labels = ['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds']
            for lab in labels:
                best_model.update({lab: float(unpack_options(cfgsetup, 'Modelling', lab)[3])})
                samples_file.update({lab: np.atleast_1d(float(unpack_options(cfgsetup, 'Modelling', lab)[3]))})

            best_model.update({'chi2': 0})
            best_model.update({'chi2/dof': 0})
            best_model.update({'id': -1})
            best_model.update({'chain': -1})
            best_model.update({'date_save': -1})
            best_model.update({'time_save': -1})
            best_model.update({'accrate': 0})
            best_model.update({'fullid': -1})

        # -------------------------------------------------------------------
        def lens_rotation(alpha0, s0, dalpha, ds, t, tb):
            '''
            Compute the angle and the primary-secondary distance with a lens
            orbital motion.

            :alpha0: angle at tb
            :s0: separation at tb
            :dalpha: lens rotation velocity in rad.year^-1
            :ds: separation velocity in year^-1
            :t: numpy array of observation dates
            :tb: a reference date
            :return: numpy array of the value of the angle and separation
                     for each date
            '''
            Cte_yr_d = 365.25  # Julian year in days
            alpha = alpha0 - (t - tb) * dalpha / Cte_yr_d  # (-1) because if the caustic rotates of dalpha, it is as if the source follows a trajectory with an angle of alpha(tb)-dalpha.
            s = s0 + (t - tb) * ds / Cte_yr_d

            return alpha, s

        # Parameters contraction
        s0 = best_model['s']
        q = best_model['q']
        u0 = best_model['u0']
        alpha0 = best_model['alpha']
        tE = best_model['tE']
        t0 = best_model['t0']
        piEN = best_model['piEN']
        piEE = best_model['piEE']
        gamma = best_model['gamma']
        dalpha = best_model['dalpha']
        ds = best_model['ds']
        GL1 = s0 * q / (1 + q)
        GL2 = - s0 / (1 + q)
        tb = cfgsetup.getfloat('Modelling', 'tb')

        chi2_flux = 0
        chi2dof_flux = 0

        #    Best model for data
        # ------------------------------------------------------------------
        # Calculation of the amplification
        param_model = best_model
        observatories = np.unique(time_serie['obs'])
        models_lib = np.unique(time_serie['model'])
        if cfgsetup.getboolean('Plotting', 'Data'):

            time_serie.update({'x': np.empty(len(time_serie['dates']))})
            time_serie.update({'y': np.empty(len(time_serie['dates']))})

            for j in xrange(len(observatories)):
                cond2 = (time_serie['obs'] == observatories[j])

                if flag_fix_gamma:
                    param_model.update({'gamma': time_serie['gamma'][cond2][0]})

                for i in xrange(models_lib.shape[0]):
                    cond = (time_serie['model'] == models_lib[i]) &\
                           (time_serie['obs'] == observatories[j])

                    if cond.sum() > 0:
                        time_serie_export = time_serie['dates'][cond]
                        DsN_export = time_serie['DsN'][cond]
                        DsE_export = time_serie['DsE'][cond]

                        Ds_export = dict({'N': DsN_export, 'E': DsE_export})

                        try:
                            kwargs_method = dict(cfgsetup.items(models_lib[i]))
                        except:
                            kwargs_method = dict()

                        amp = models[models_lib[i]].magnifcalc(time_serie_export,
                                param_model, Ds=Ds_export, tb=tb, **kwargs_method)

                        time_serie['amp'][cond] = amp
                        # print amp
                        # print time_serie['amp'][cond]
                        del amp

                # Calculation of fs and fb
                fs, fb = fsfb(time_serie, cond2, blending=True)
                # if (fb/fs < 0 and observatories[j]=="ogle-i"):
                #     fs, fb = fsfb(time_serie, cond2, blending=False)

                time_serie['fs'][cond2] = fs
                time_serie['fb'][cond2] = fb
                # print fs, fb

                if (observatories[j] == cfgsetup.get('Observatories', 'Reference').lower()) \
                        | (j == 0):
                    fs_ref = fs
                    fb_ref = fb
                    mag_baseline = 18.0 - 2.5 * np.log10(1.0 * fs_ref + fb_ref)
                    # print fs, fb

                # Source postion
                if cond2.sum() > 0:
                    DsN = time_serie['DsN'][cond2]
                    DsE = time_serie['DsE'][cond2]

                    t = time_serie['dates'][cond2]

                    tau = (t - t0) / tE + piEN * DsN + piEE * DsE
                    beta = u0 + piEN * DsE - piEE * DsN

                    z = (tau + 1j * beta) * np.exp(1j * alpha0)

                    time_serie['x'][cond2] = z.real
                    time_serie['y'][cond2] = z.imag

            # Calculation of chi2
            time_serie['flux_model'] = time_serie['amp'] * time_serie['fs'] + \
                                       time_serie['fb']

            # time_serie['chi2pp'] = np.power((time_serie['flux'] - time_serie[
            #     'flux_model']) / time_serie['err_flux'], 2)
            time_serie['residus'] = time_serie['magnitude'] - (18.0 - 2.5 * np.log10(time_serie['flux_model']))
            time_serie['residus_flux'] = time_serie['flux'] - time_serie['flux_model']
            time_serie['mgf_data'] = (time_serie['flux'] - time_serie['fb']) / time_serie['fs']
            time_serie['mgf_data_err'] = time_serie['err_flux'] / time_serie['fs']
            time_serie['res_mgf'] = time_serie['mgf_data'] - time_serie['amp']
            time_serie['chi2pp'] = np.power(time_serie['residus'] / time_serie['err_magn'], 2.0)
            time_serie['chi2pp_flux'] = np.power(time_serie['residus_flux'] / time_serie['err_flux'], 2.0)
            chi2 = np.sum(time_serie['chi2pp'])
            chi2_flux = np.sum(time_serie['chi2pp_flux'])

        # Calculation of the lightcurve model
        plot_min = float(options.split('/')[0].split('-')[0].strip())
        plot_max = float(options.split('/')[0].split('-')[1].strip())
        nb_pts = int(options.split('/')[1].strip())

        locations = np.unique(obs_properties['loc'])


        print "      Max magnification: {:.2f}".format(np.max(time_serie['amp']))
        # print len(time_serie['amp'])

        # Fit summary
        text = "Fit summary"
        communicate(cfgsetup, 3, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)
        observatories_com = np.unique(time_serie['obs'])

        if (cfgsetup.getint("Modelling", "Verbose") >= 3) & (rang_best_model == rang_2plot[idmod]):
            fn_output_terminal = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Outputs')\
                                 + "Results.txt"
            # print fn_output_terminal
            file = open(fn_output_terminal, 'w')

            text = "\n\033[1m\033[7m  {:25s} {:>9s}    {:>9s}    {:>9s}    {:>9s}    \033[0m".format(
                "Site", "chi^2", "chi^2/dof", "RF 1", "RF 2")
            print text

            text_precis = "\n  {:25s} {:>18s}    {:>18s}    {:>18s}    {:>18s}    \n".format(
                "Site", "chi^2", "chi^2/dof", "RF 1", "RF 2")
            file.write(text_precis)

            params_raw = np.array(['t0', 'u0', 'tE', 'rho', 'gamma', 'piEN', 'piEE', 's', 'q', 'alpha', 'dalpha', 'ds'])
            n_param = nb_param_fit

            # n_param = int(len(time_serie['dates']) - chi2 / samples_file['chi2/dof'][rang_best_model])
            nb_data_tot = 0
            # observatories_com = np.unique(time_serie['obs'])
            text = ""
            text2 = ""
            text_precis = ""
            text2_precis = ""
            for i in xrange(len(observatories_com)):
                rf = [float(a.replace("(", "").replace(")", "").strip()) for a in unpack_options(cfgsetup, "Observatories", observatories_com[i])[:2]]
                cond = time_serie['obs']==observatories_com[i]
                chi2_com = np.sum(time_serie['chi2pp'][cond])
                nb_data_tot = nb_data_tot + len(time_serie['chi2pp'][cond])
                if len(time_serie['chi2pp'][cond]) > 0:
                    chi2dof_com = chi2_com / (len(time_serie['chi2pp'][cond])-n_param)
                else:
                    chi2dof_com = 0.0
                text = text + "  {:25s} {:9.3e}    {:9.3e}    {:9.3e}    {:9.3e}\n".format(observatories_com[i].upper(), chi2_com, chi2dof_com, rf[0], rf[1])
                text_precis = text_precis + "  {:25s} {:18.12e}    {:18.12e}    {:18.12e}    {:18.12e}\n".format(observatories_com[i].upper(), chi2_com, chi2dof_com, rf[0], rf[1])

                try:
                    g = time_serie["fb"][cond][0]/time_serie["fs"][cond][0]
                except:
                    g = np.inf
                # Y = 18.0 - 2.5 * np.log10(time_serie['fb'][cond][0] + time_serie['fs'][cond][0])
                Y = time_serie['fb'][cond][0] + time_serie['fs'][cond][0]
                # Yb = 18.0 - 2.5 * np.log10(time_serie['fb'][cond][0])
                Yb = time_serie['fb'][cond][0]
                # Ys = 18.0 - 2.5 * np.log10(time_serie['fs'][cond][0])
                Ys = time_serie['fs'][cond][0]
                text2 = text2 + "  {:25s} {:8.3f}   {:8.3f}   {:8.3f}   {:8.3f}   {:5.3f}\n".format(
                    observatories_com[i].upper(), Y, Yb, Ys, g, gamma)

                text2_precis = text2_precis + "  {:25s} {:18.12e}   {:18.12e}   {:18.12e}   {:18.12e}   {:18.12e}\n".format(
                    observatories_com[i].upper(), Y, Yb, Ys, g, gamma)

                if (observatories[i] == cfgsetup.get('Observatories', 'Reference').lower()) \
                        | (i == 0):

                    Y = 18.0 - 2.5 * np.log10(time_serie['fb'][cond][0] + time_serie['fs'][cond][0])
                    Yb = 18.0 - 2.5 * np.log10(time_serie['fb'][cond][0])
                    Ys = 18.0 - 2.5 * np.log10(time_serie['fs'][cond][0])
                    text3 = "Reference for magnitudes:\n  {:25s} {:8.3f}   {:8.3f}   {:8.3f}\n".format(
                            observatories_com[i].upper(), Y, Yb, Ys)

                    text3_precis = "Reference for magnitudes:\n  {:25s} {:18.12e}   {:18.12e}   {:18.12e}\n".format(
                        observatories_com[i].upper(), Y, Yb, Ys)

            print text
            file.write(text_precis)

            text = "{:25}={:2}{:9.3e}{:4}{:9.3e} (chi^2 on magn)".format("", "", chi2_flux, "", chi2/(nb_data_tot-n_param))
            chi2dof_flux = chi2/(nb_data_tot-n_param)
            print text
            text = "{:25}={:2}{:18.12e}{:4}{:18.12e} (chi^2 on magn)".format("", "", chi2_flux, "", chi2/(nb_data_tot-n_param))
            chi2dof_flux = chi2/(nb_data_tot-n_param)
            file.write(text)

            text = "\n\033[1m\033[7m  {:78s}\033[0m".format("Best-fitting parameters")
            print text
            text = "\n  {:78s}\n".format("Best-fitting parameters")
            file.write(text)

            #print samples_file
            piE = np.sqrt(np.power(samples_file['piEN'][rang_best_model], 2) + np.power(samples_file['piEE'][rang_best_model],2))
            gamma = np.sqrt((samples_file['ds'][rang_best_model]/samples_file['s'][rang_best_model])**2 + samples_file['dalpha'][rang_best_model]**2)
            text = "{:>10} = {:.6f}\n".format("q", samples_file['q'][rang_best_model]) + "{:>10} = {:.6f}\n".format("s", samples_file['s'][rang_best_model]) + "{:>10} = {:.6f}\n".format("tE", samples_file['tE'][rang_best_model]) + "{:>10} = {:.6f}\n".format("rho", samples_file['rho'][rang_best_model]) + "{:>10} = {:.6f}\n".format("piEN", samples_file['piEN'][rang_best_model]) + "{:>10} = {:.6f}\n".format("piEE", samples_file['piEE'][rang_best_model]) + "{:>10} = {:.6f}\n".format("piE", piE) + "{:>10} = {:.6f}\n".format("t0", samples_file['t0'][rang_best_model]) + "{:>10} = {:.6f}\n".format("u0", samples_file['u0'][rang_best_model]) + "{:>10} = {:.6f}\n".format("alpha", samples_file['alpha'][rang_best_model]) + "{:>10} = {:.6f}\n".format("dalpha", samples_file['dalpha'][rang_best_model]) + "{:>10} = {:.6f}\n".format("ds", samples_file['ds'][rang_best_model]) + "{:>10} = {:.6f}\n".format("gammaL", gamma) + "{:>10} = {:.6f}\n".format("tp", cfgsetup.getfloat("Modelling", "tp")) + "{:>10} = {:.6f}\n".format("tb", cfgsetup.getfloat("Modelling", "tb"))
            print text
            text = "{:>10} = {:.12e}\n".format("q", samples_file['q'][rang_best_model]) + "{:>10} = {:.12e}\n".format("s", samples_file['s'][rang_best_model]) + "{:>10} = {:.12e}\n".format("tE", samples_file['tE'][rang_best_model]) + "{:>10} = {:.12e}\n".format("rho", samples_file['rho'][rang_best_model]) + "{:>10} = {:.12e}\n".format("piEN", samples_file['piEN'][rang_best_model]) + "{:>10} = {:.12e}\n".format("piEE", samples_file['piEE'][rang_best_model]) + "{:>10} = {:.12e}\n".format("piE", piE) + "{:>10} = {:.12e}\n".format("t0", samples_file['t0'][rang_best_model]) + "{:>10} = {:.12e}\n".format("u0", samples_file['u0'][rang_best_model]) + "{:>10} = {:.12e}\n".format("alpha", samples_file['alpha'][rang_best_model]) + "{:>10} = {:.12e}\n".format("dalpha", samples_file['dalpha'][rang_best_model]) + "{:>10} = {:.12e}\n".format("ds", samples_file['ds'][rang_best_model]) + "{:>10} = {:.12e}\n".format("gammaL", gamma) + "{:>10} = {:.12e}\n".format("tp", cfgsetup.getfloat("Modelling", "tp")) + "{:>10} = {:.12e}\n".format("tb", cfgsetup.getfloat("Modelling", "tb"))
            file.write(text)

            text = "\n\033[1m\033[7m  {:25s} {:>8s}   {:>8s}   {:>8s}     {:>6s}   {:>5s}{:3s}\033[0m".format(
                "Site", "Baseline", "Blending", "Source", "Fb/Fs", "LLD", "")
            print text
            text = "\n  {:25s} {:>18s}   {:>18s}   {:>18s}   {:>18s}   {:>18s}{:3s}\n".format(
                "Site", "Baseline", "Blending", "Source", "Fb/Fs", "LLD", "")
            file.write(text)
            print text2
            file.write(text2_precis)
            print text3
            file.write(text3_precis)
            file.close()

        #    Best model theoretical light curve
        # ------------------------------------------------------------------
        model2load = np.array([])
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
            'RelativePaths', 'Data')

        if len(obs_properties['loc']) > 1:
            name1 = obs_properties['loc'][np.where(np.array(
                [obs == cfgsetup.get('Observatories', 'Reference').lower()
                 for obs in observatories]) == True)[0][0]]
            name1 = glob.glob(path + name1 + '.*')[0]
        else:
            name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]

        for i in xrange(len(locations)):
            name = 'Models_' + locations[i]
            models_temp = model_param[name]
            name = 'DateRanges_' + locations[i]
            dates_temp = model_param[name]

            # min = []
            # max = []
            # for j in xrange(len(models_temp)):
            #       model2load = np.append(model2load, models_temp[j])
            #       tmin = float((dates_temp[j]).split('-')[0].strip())
            #       tmax = float((dates_temp[j]).split('-')[1].strip())
            #
            #       min = np.append(min, tmin)
            #       max = np.append(min, tmax)
            #
            # min = np.min(min)
            # max = np.max(max)

            min = plot_min
            max = plot_max

            if i == 0:
                model_time_serie = np.array([dict({
                    'dates': np.linspace(min, max, nb_pts),
                    'model': np.full(nb_pts, '0', dtype='S100'),
                    'amp': np.full(nb_pts, 0.1, dtype='f8'),
                })])
            else:
                model_time_serie = np.append(model_time_serie, np.array([dict({
                    'dates': np.linspace(min, max, nb_pts),
                    'model': np.full(nb_pts, '0', dtype='S100'),
                    'amp': np.full(nb_pts, 0.1, dtype='f8'),
                })]))

            for j in xrange(len(models_temp)):
                tmin = float((dates_temp[j]).split('-')[0].strip())
                tmax = float((dates_temp[j]).split('-')[1].strip())

                cond = (model_time_serie[i]['dates'] <= tmax) \
                       & (model_time_serie[i]['dates'] >= tmin)

                model_time_serie[i]['model'][cond] = models_temp[j]

            cond = model_time_serie[i]['model'] == '0'
            if cond.sum() > 0:
                model_time_serie[i]['model'][cond] = models_temp[0]

            # Ephemeris
            c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), \
                              dec=cfgsetup.get('EventDescription', 'DEC'),
                              frame='icrs')
            # print c_icrs.transform_to('barycentrictrueecliptic')
            l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
            b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

            name2 = glob.glob(path + locations[i] + '.*')[0]
            sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
                ephemeris.Ds(name1, name2, l, b,
                             cfgsetup.getfloat('Modelling', 'tp'), \
                             cfgsetup)

            if name1 != name2:
                DsN = DsNs
                DsE = DsEs
            else:
                DsN = DsNe
                DsE = DsEe

            model_time_serie[i].update({'DsN': np.array(
                [DsN(a) for a in model_time_serie[i]['dates']])})
            model_time_serie[i].update({'DsE': np.array(
                [DsE(a) for a in model_time_serie[i]['dates']])})

            # Amplification
            models_lib = np.unique(model_time_serie[i]['model'])
            for k in xrange(models_lib.shape[0]):
                cond = (model_time_serie[i]['model'] == models_lib[k])

                if cond.sum() > 0:
                    time_serie_export = model_time_serie[i]['dates'][cond]
                    DsN_export = model_time_serie[i]['DsN'][cond]
                    DsE_export = model_time_serie[i]['DsE'][cond]

                    Ds_export = dict({'N': DsN_export, 'E': DsE_export})

                    try:
                        kwargs_method = dict(cfgsetup.items(models_lib[k]))
                    except:
                        kwargs_method = dict()

                    amp = models[models_lib[k]].magnifcalc(time_serie_export, param_model, Ds=Ds_export, tb=tb, **kwargs_method)

                    model_time_serie[i]['amp'][cond] = amp

            if cfgsetup.getboolean('Plotting', 'Data'):
                model_time_serie[i].update({'magnitude': 18.0 - 2.5 * np.log10(
                    fs_ref * model_time_serie[i]['amp'] + fb_ref)})

            # Source position in (x, y)
            DsN = model_time_serie[i]['DsN']
            DsE = model_time_serie[i]['DsE']

            t = model_time_serie[i]['dates']

            tau = (t - t0) / tE + piEN * DsN + piEE * DsE
            beta = u0 + piEN * DsE - piEE * DsN

            z = (tau + 1j * beta) * np.exp(1j * alpha0)

            model_time_serie[i].update({'x': z.real})
            model_time_serie[i].update({'y': z.imag})

        del amp, DsN_export, DsE_export, Ds_export, cond, time_serie_export

        # print model_time_serie[1]['model']

        # # Interpolation method
        # # -------------------------------------------------------------------------
        # key_list = [key for key in interpol_method]
        #
        # interpol_func = dict()
        # if len(key_list) > 0:
        #     for i in xrange(len(key_list)):
        #         time_serie_export = interpol_method[key_list[i]][0]
        #
        #         DsN_export = interpol_method[key_list[i]][1]
        #         DsE_export = interpol_method[key_list[i]][2]
        #
        #         Ds_export = dict({'N':DsN_export, 'E':DsE_export})
        #
        #         name = key_list[i].split('#')[1]
        #         amp = models[name].magnifcalc(time_serie_export, param_model, Ds=Ds_export)
        #
        #         interpol_method[key_list[i]][3] = amp
        #         # interpol_func = interpolate.interp1d(time_serie_export, amp)
        #         interpol_func.update({key_list[i]: interpolate.interp1d(time_serie_export, amp, kind='linear')})

        #    Reference frames
        # ------------------------------------------------------------------
        # Orientations on the Sky
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
            'RelativePaths', 'Data')
        if len(model_time_serie) == 2:
            EarthSun, EarthSat, D_enm = boussole(
                EarthSunFile=path + "Earth.dat",
                EarthSatelliteFile=path + "Spitzer.dat",
                cfg=cfgsetup, t_D_xy=best_model['t0'])
        else:
            EarthSun, EarthSat, D_enm = boussole(
                EarthSunFile=path + "Earth.dat",
                EarthSatelliteFile=path + "Earth.dat",
                cfg=cfgsetup, t_D_xy=best_model['t0'])

        #    Sigma clipping
        # ------------------------------------------------------------------
        # Determine the best rescaling factors
        if (cfgsetup.getint("Modelling", "Verbose") > 4) & (nb_param_fit > 0):
            text = "\n\033[1m\033[7m{:>2s}{:<25s}{:1s}{:>10s}{:1s}{:>5s}{:1s}{:>10s}{:1s}{:>5s}{:1s}{:>10s}{:1s}{:>5s}{:2s}\033[0m".format(
                    "", "Site", "", "RF1(loop3)", "", "Rej.", "", "RF1(loop5)", "", "Rej.", "", "RF1(loop7)", "", "Rej.", "")
            print text

            def func(f1, table, f2, ddl):
                x = np.sum(np.power(table['residus'], 2)/(np.power(f1*table['err_magn'], 2) + f2**2))
                x = x / ddl - 1.0
                return x


            text = ""
            for j in xrange(len(observatories_com)):
                # Pre-defied rescaling factors
                f1 = float(unpack_options(cfgsetup, 'Observatories', observatories[0])[0].replace('(', ''))
                f2 = float(unpack_options(cfgsetup, 'Observatories', observatories[0])[1].replace(')', ''))
                if abs(f1-1.0) > 1e-10:
                    text = "{:>2s}{:<25s}{:<30s}\n".format("", observatories_com[j].upper(), "RF 1 not equal to 1.0.")
                    continue
                # Select the observatory
                condj = np.where(time_serie['obs'] == observatories[j])
                time_serie_SC = copy.deepcopy(time_serie)
                [time_serie_SC.update({key: time_serie_SC[key][condj]}) for key in time_serie_SC]
                # Compute the degree of freedom ddl
                nb_data = len(time_serie_SC['dates'])
                if nb_data > nb_param_fit:
                    ddl = nb_data - nb_param_fit
                else:
                    ddl = nb_data

                # Compute the rescaling factor f1 from the value of f2
                rejected_points_id = np.array([])
                nb_reject_sc = 0
                nb_loops = 7
                text = text + "{:>2s}{:<25s}".format("", observatories_com[j].upper())
                for i in xrange(nb_loops):
                    mean = np.mean(time_serie_SC['err_magn'])
                    sdt = np.std(time_serie_SC['err_magn'])
                    toremove = np.where(np.abs(time_serie_SC['err_magn'] - mean) > 3.0 * sdt)
                    nb_reject_sc = nb_reject_sc + len(toremove[0])
                    if len(toremove[0]) > 0:
                        rejected_points_id = np.append(rejected_points_id, time_serie_SC['id'][toremove])
                        [time_serie_SC.update({key : np.delete(time_serie_SC[key], toremove)}) for key in time_serie_SC]

                    if (i==2) | (i==4) | (i==6):
                        try:
                            f1_op = fsolve(func, 1.0, args=(time_serie_SC, f2, ddl))
                        except:
                            f1_op = 0.0
                        text = text + "{:>10.3f}{:1s}{:>5d}{:1s}".format(
                                f1_op[0], "", nb_reject_sc, "")

                text = text + "\n"
            print text


        # ---------------------------------------------------------------------
        #    Create an html webpage (amplification if no data, magnitude if so)
        # ---------------------------------------------------------------------
        if cfgsetup.getboolean('Plotting', 'Data'):
            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                'RelativePaths', 'Plots')

            if (best_model['fullid'] == -1) & flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-summary.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ': best model last MCMC'
            elif flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-summary-' \
                           + repr(best_model['fullid']) + '.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model # ' + repr(
                    best_model['fullid'])
            else:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-summary-fix.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model fix'

            if os.path.exists(filename):
                os.remove(filename)
            bplt.output_file(filename)
            fig = np.array([])

            # Preparation of the data
            time_serie.update({'colour': np.full(len(time_serie['dates']), 'black',
                                                 dtype='S100')})
            time_serie.update(
                {'mag_align': np.full(len(time_serie['dates']), 0, dtype='f8')})

            palette = plt.get_cmap('Blues')

            observatories = np.unique(time_serie['obs'])
            for i in xrange(len(observatories)):
                cond = np.where(time_serie['obs'] == observatories[i])

                cond2 = \
                np.where(observatories[i] == np.array(obs_properties['key']))[0][0]
                color = '#' + obs_properties['colour'][cond2]

                time_serie['colour'][cond] = color

                # Magnitude aligned
                Y = ((time_serie['flux'][cond] - time_serie['fb'][cond])/time_serie['fs'][cond])    #############
                Y = 18.0 - 2.5 * np.log10(fs_ref * Y + fb_ref)
                # cond3 = cond & (Y>0)                                                           ##################
                # Y[cond3] = 18.0 - 2.5 * np.log10(fs_ref * Y[cond3] + fb_ref)                 ##################
                # cond3 = cond & (Y<=0)                                               ##################
                # Y[cond3] = 1000                                                    ##################
                time_serie['mag_align'][cond] = Y                               ##################


            #    Create output files
            # ---------------------------------------------------------------------
            path_outputs = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Outputs')
            if not os.path.exists(path_outputs):
                os.makedirs(path_outputs)

            for j in xrange(len(observatories_com)):
                
                idx = [jj for jj in xrange(len(observatories_com)) if observatories_com[j]==obs_properties['key'][jj]][0]
                flag_fom = obs_properties['fluxoumag'][idx]

                if flag_fom.lower()=='magnitude':
                    text = "#{11:>5s} {0:>18s} {1:>6s} {3:>12s} {4:>10s} {8:>8s} {9:>9s} {10:>9s} {5:>12s} {6:>12s} {13:>10s} {14:>10s} {7:>9s} {12:>20s} {2:>20s}\n".format(
                            "Date", "Magn", "Err_Magn", "Err_Magn_Res", "Resi", "Back", "Seeing", "Chi2", "Mgf-dat", "Err_Mgf", "Resi-Mgf", "ID", "Input_Magn", "x", "y")

                elif flag_fom.lower()=='flux':
                    text = "#{11:>5s} {0:>18s} {1:>6s} {3:>12s} {4:>10s} {8:>8s} {9:>9s} {10:>9s} {5:>12s} {6:>12s} {13:>10s} {14:>10s} {7:>9s} {12:>20s} {2:>20s}\n".format(
                            "Date", "Magn", "Err_Flux", "Err_Magn_Res", "Resi", "Back", "Seeing", "Chi2", "Mgf-dat", "Err_Mgf", "Resi-Mgf", "ID", "Input_Flux")

                filename = path_outputs + observatories_com[j].upper() + ".dat"

                condj = np.where(time_serie['obs'] == observatories[j])
                time_serie_SC = copy.deepcopy(time_serie)
                [time_serie_SC.update({key: time_serie_SC[key][condj]}) for key in time_serie_SC]

                if flag_fom.lower()=='magnitude':
                    for jj in xrange(len(time_serie_SC['dates'])):
                        text = text +\
                                "{11:6d} {0:18.12f} {1:6.3f} {3:12.3e} {4:10.3e} {8:8.3f} {9:9.3e} {10:9.2e} {5:12.5f} {6:12.5f} {13:10.6f} {14:10.6f} {7:9.3e} {12:20.12f} {2:20.12f}".format(
                                time_serie_SC['dates'][jj],
                                time_serie_SC['mag_align'][jj],
                                time_serie_SC['err_magn_orig'][jj],
                                time_serie_SC['err_magn'][jj],
                                time_serie_SC['residus'][jj],
                                time_serie_SC['background'][jj],
                                time_serie_SC['seeing'][jj],
                                time_serie_SC['chi2pp'][jj],
                                time_serie_SC['mgf_data'][jj],
                                time_serie_SC['mgf_data_err'][jj],
                                time_serie_SC['res_mgf'][jj],
                                time_serie_SC['id'][jj],
                                time_serie_SC['magnitude'][jj],
                                time_serie_SC['x'][jj],
                                time_serie_SC['y'][jj]
                                )
                        text = text + "\n"

                elif flag_fom.lower()=='flux':
                    for jj in xrange(len(time_serie_SC['dates'])):
                        text = text +\
                                "{11:6d} {0:18.12f} {1:6.3f} {3:12.3e} {4:10.3e} {8:8.3f} {9:9.3e} {10:9.2e} {5:12.5f} {6:12.5f} {13:10.6f} {14:10.6f} {7:9.3e} {12:20.12f} {2:20.12f}".format(
                                time_serie_SC['dates'][jj],
                                time_serie_SC['mag_align'][jj],
                                time_serie_SC['err_flux_orig'][jj],
                                time_serie_SC['err_magn'][jj],
                                time_serie_SC['residus'][jj],
                                time_serie_SC['background'][jj],
                                time_serie_SC['seeing'][jj],
                                time_serie_SC['chi2pp'][jj],
                                time_serie_SC['mgf_data'][jj],
                                time_serie_SC['mgf_data_err'][jj],
                                time_serie_SC['res_mgf'][jj],
                                time_serie_SC['id'][jj],
                                time_serie_SC['flux'][jj],
                                time_serie_SC['x'][jj],
                                time_serie_SC['y'][jj]
                                )
                        text = text + "\n"

                file = open(filename, 'w')
                file.write(text)
                file.close()




            # ..................................................................
            #    Plot light curve : plc
            # ..................................................................

            cond = np.isnan(time_serie['mag_align'])
            time_serie['mag_align'][cond] = 0

            source = ColumnDataSource(time_serie)
            col_used = ['id', 'obs', 'dates', 'mag_align', 'x', 'y', 'colour', 'residus']
            col_all = [key for key in time_serie]
            [col_all.remove(key) for key in col_used]
            [source.remove(key) for key in col_all]
            hover_plc = HoverTool(tooltips=[("ID", "@id{int}"), ("Obs", "@obs"), ("Date", "@dates{1.11}")])

            tmin = float(options.split('/')[0].split('-')[0].strip())
            tmax = float(options.split('/')[0].split('-')[1].strip())

            ymin = np.min(time_serie['mag_align'])
            ymax = np.max(time_serie['mag_align'])

            tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap", hover_plc]
            fig = np.append(fig, \
                            bplt.figure(toolbar_location="above", plot_width=1200,
                                        plot_height=600, x_range=(tmin, tmax),
                                        y_range=(ymax, ymin), \
                                        title=None, min_border=10,
                                        min_border_left=50, tools=tools))

            fig_curr = fig[0]

            # Annotations
            # ^^^^^^^^^^^
            colours = ['black', '#297CC4']
            id_colour = 0
            for i in xrange(len(locations)):
                name = 'Models_' + locations[i]
                models_temp = model_param[name]
                name = 'DateRanges_' + locations[i]
                dates_temp = model_param[name]

                # print np.max(model_time_serie[i]['amp'])

                for j in xrange(len(models_temp)):
                    tmin = float((dates_temp[j]).split('-')[0].strip())
                    tmax = float((dates_temp[j]).split('-')[1].strip())

                    X = np.ones(2) * tmin
                    Y = np.linspace(0, 100000, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                    X = np.ones(2) * tmax
                    Y = np.linspace(0, 100000, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Amplification models
            # ^^^^^^^^^^^^^^^^^^^^
            if cfgsetup.getboolean("Plotting", "Data"):
                colours = ['black', '#297CC4']
                id_colour = 0
                for i in xrange(len(locations)):
                    X = model_time_serie[i]['dates']
                    Y = model_time_serie[i]['magnitude']
                    fig_curr.line(X, Y, line_width=2, color=colours[id_colour],
                                  alpha=1)

                    if id_colour < len(colours) - 1:
                        id_colour = id_colour + 1
                    else:
                        id_colour = 0

                # Write output files for the models
                text = "#{0:>17s} {1:>9s} {2:>6s} {3:>7s} {4:>7s}\n".format("Date", "Mgf", "Magn", "x", "y")
                filename = path_outputs + locations[i].upper() + ".dat"

                time_serie_SC = copy.deepcopy(model_time_serie[i])
                [time_serie_SC.update({key: time_serie_SC[key]}) for key in time_serie_SC]

                for jj in xrange(len(time_serie_SC['dates'])):
                    text = text +\
                            "{0:18.12f} {1:9.3f} {2:6.3f} {3:7.3f} {4:7.3f}".format(
                            time_serie_SC['dates'][jj],
                            time_serie_SC['amp'][jj],
                            time_serie_SC['magnitude'][jj],
                            time_serie_SC['x'][jj],
                            time_serie_SC['y'][jj]
                            )
                    text = text + "\n"

                file = open(filename, 'w')
                file.write(text)
                file.close()

                # Magnitude
                fig_curr.circle('dates', 'mag_align', size=8, color='colour',
                                alpha=0.4, source=source)


                # Legend
                for i in xrange(len(obs_properties['name'])):
                    col = '#' + obs_properties['colour'][i]
                    fig_curr.circle(-10000, -10000, size=8, color=col, alpha=0.4,
                                    legend=obs_properties['name'][i])

                # Magnitude (errors)
                # ^^^^^^^^^^^^^^^^^^
                err_xs = []
                err_ys = []

                for x, y, yerr, colori in zip(time_serie['dates'],
                                              time_serie['mag_align'],
                                              time_serie['err_magn'],
                                              time_serie['colour']):
                    err_xs.append((x, x))
                    err_ys.append((y - yerr, y + yerr))

                fig_curr.multi_line(err_xs, err_ys, color=time_serie['colour'])

            # Layout
            # ^^^^^^
            fig_curr.xaxis.axis_label = 'HJD - 2,450,000'
            fig_curr.yaxis.axis_label = 'Magnitude'
            fig_curr.xaxis.axis_label_text_font = 'helvetica'
            fig_curr.yaxis.axis_label_text_font = 'helvetica'
            fig_curr.xaxis.axis_label_text_font_size = '10pt'
            fig_curr.yaxis.axis_label_text_font_size = '10pt'

            fig_curr.min_border_top = 10
            fig_curr.min_border_bottom = 0
            fig_curr.min_border_left = 0

            fig_curr.xgrid.grid_line_color = None
            fig_curr.ygrid.grid_line_color = None

            # ..................................................................
            #    Plot residus in mag : prm
            # ..................................................................

            hover_prm = HoverTool(
                tooltips=[
                    ("ID", "@id{int}"),
                    ("Obs", "@obs"),
                    ("Date", "@dates{1.11}")
                ]
            )

            tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap",
                     hover_prm]
            fig = np.append(fig, \
                            bplt.figure(toolbar_location="above", plot_width=1200,
                                        plot_height=300, x_range=fig[0].x_range,
                                        y_range=(-0.25, 0.25), \
                                        title=None, min_border=10,
                                        min_border_left=50, tools=tools))

            fig_curr = fig[1]

            # Annotations
            # ^^^^^^^^^^^
            colours = ['black', '#297CC4']
            id_colour = 0
            for i in xrange(len(locations)):
                name = 'Models_' + locations[i]
                models_temp = model_param[name]
                name = 'DateRanges_' + locations[i]
                dates_temp = model_param[name]

                for j in xrange(len(models_temp)):
                    tmin = float((dates_temp[j]).split('-')[0].strip())
                    tmax = float((dates_temp[j]).split('-')[1].strip())

                    X = np.ones(2) * tmin
                    Y = np.linspace(-100, 100, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                    X = np.ones(2) * tmax
                    Y = np.linspace(-100, 100, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Magnitude
            fig_curr.circle('dates', 'residus', size=8, color='colour', alpha=0.4,
                            source=source)

            # Magnitude (errors)
            # ^^^^^^^^^^^^^^^^^^
            err_xs = []
            err_ys = []

            for x, y, yerr, colori in zip(time_serie['dates'],
                                          time_serie['residus'],
                                          time_serie['err_magn'],
                                          time_serie['colour']):
                err_xs.append((x, x))
                err_ys.append((y - yerr, y + yerr))

            fig_curr.multi_line(err_xs, err_ys, color=time_serie['colour'])

            X = np.linspace(-100000, 100000, 2)
            Y = np.ones(2) * 0
            fig_curr.line(X, Y, line_width=1, color='dimgray', alpha=1)

            X = np.linspace(-100000, 100000, 2)
            Y = np.ones(2) * 0.1
            fig_curr.line(X, Y, line_width=0.5, color='dimgray', alpha=0.5)

            X = np.linspace(-100000, 100000, 2)
            Y = np.ones(2) * (-0.1)
            fig_curr.line(X, Y, line_width=0.5, color='dimgray', alpha=0.5)

            # Layout
            # ^^^^^^
            fig_curr.xaxis.axis_label = 'HJD - 2,450,000'
            fig_curr.yaxis.axis_label = 'Residuals [mag]'
            fig_curr.xaxis.axis_label_text_font = 'helvetica'
            fig_curr.yaxis.axis_label_text_font = 'helvetica'
            fig_curr.xaxis.axis_label_text_font_size = '10pt'
            fig_curr.yaxis.axis_label_text_font_size = '10pt'

            fig_curr.min_border_top = 10
            fig_curr.min_border_bottom = 0
            fig_curr.min_border_left = 0

            fig_curr.xgrid.grid_line_color = None
            fig_curr.ygrid.grid_line_color = None

            # ..................................................................
            #    Plot caustic 1 : pc1
            # ..................................................................

            hover_pc1 = HoverTool(
                tooltips=[
                    ("ID", "@id{int}"),
                    ("Obs", "@obs"),
                    ("Date", "@dates{1.11}")
                ]
            )

            xmin = -1.0
            xmax = 1.0
            ymin = -1.0
            ymax = 1.0
            tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap",
                     hover_pc1]
            fig = np.append(fig, \
                            bplt.figure(toolbar_location="above", plot_width=600,
                                        plot_height=560, x_range=(xmin, xmax),
                                        y_range=(ymin, ymax), \
                                        title=None, min_border=10,
                                        min_border_left=50, tools=tools))

            fig_curr = fig[2]

            # Caustic
            # ^^^^^^^

            # Case of lens orbital rotation
            try:
                time_caustic = options.split('/')[2].replace('[','').replace(']','').split('-')
                time_caustic = np.array([float(a.strip()) for a in time_caustic])
                n_caustics = len(time_caustic)
                nb_pts_caus = 1000
                alpha, s = lens_rotation(alpha0, s0, dalpha, ds, time_caustic, tb)
                color_caustics = np.array(['Orange', 'SeaGreen', 'LightSeaGreen', 'CornflowerBlue', 'DarkViolet'])
            except:
                nb_pts_caus = 1000
                n_caustics = 0

            if q > 1e-9:
                if n_caustics > 0:
                    numerous_caustics = []
                    for i in xrange(n_caustics):
                        # > Courbes critiques et caustiques
                        delta = 2.0 * np.pi / (nb_pts_caus - 1.0)
                        phi = np.arange(nb_pts_caus) * delta
                        critic = np.array([])
                        for phi_temp in phi[:]:
                            critic = np.append(critic, critic_roots(s[i], q, phi_temp))
                        caustic = critic - 1 / (1 + q) * (1 / critic.conjugate() + q / (critic.conjugate() + s[i]))
                        caustic = caustic + GL1  # From Cassan (2008) to CM
                        caustic = caustic * np.exp(1j*(alpha[i]-alpha0))

                        fig_curr.circle(caustic.real, caustic.imag, size=0.5, color=color_caustics[0], alpha=0.5)
                        print color_caustics
                        color_caustics = np.roll(color_caustics, -1)
                        
                        numerous_caustics.append(caustic)

                # > Courbes critiques et caustiques
                delta = 2.0 * np.pi / (nb_pts_caus - 1.0)
                phi = np.arange(nb_pts_caus) * delta
                critic = np.array([])
                for phi_temp in phi[:]:
                    critic = np.append(critic, critic_roots(s0, q, phi_temp))
                caustic = critic - (1.0/(1 + q)) * (1/critic.conjugate() + q/(critic.conjugate() + s0))
                caustic = caustic + GL1  # From Cassan (2008) to CM

                fig_curr.circle(caustic.real, caustic.imag, size=0.5, color='red',
                                alpha=0.5)
                fig_curr.circle(GL1.real, GL1.imag, size=5, color='orange', alpha=1)
                fig_curr.circle(GL2.real, GL2.imag, size=5, color='orange', alpha=1)

                # Write output files
                text = "#{:>19s} {:>20s}".format("x", "y")
                if n_caustics > 0:
                    for i in xrange(n_caustics):
                            text = text +\
                                    " x({0:17.6f}) y({0:17.6f})".format(time_caustic[i])
                text = text + "\n"

                filename = path_outputs + "CAUSTIC.dat"

                for jj in xrange(len(caustic.real)):
                    text = text +\
                            "{:20.12f} {:20.12f}".format(
                            caustic.real[jj],
                            caustic.imag[jj]
                            )

                    if n_caustics > 0:
                        for i in xrange(n_caustics):
                                text = text +\
                                        " {:20.12f} {:20.12f}".format(
                                        numerous_caustics[i].real[jj],
                                        numerous_caustics[i].imag[jj]
                                        )
                    text = text + "\n"

                file = open(filename, 'w')
                file.write(text)
                file.close()
            else:
                fig_curr.circle(0, 0, size=5, color='orange', alpha=1)

            # Trajectories
            # ^^^^^^^^^^^^
            colours = ['black', '#297CC4']
            id_colour = 0

            for i in xrange(len(locations)):
                temp = np.array([abs(a - best_model['t0']) for a in
                                 model_time_serie[i]['dates']])
                rang_c = np.where(temp == np.min(temp))[0][0]

                X = model_time_serie[i]['x']
                Y = model_time_serie[i]['y']
                fig_curr.line(X, Y, line_width=2, color=colours[id_colour],
                              alpha=1)

                # Arrows
                n = 0
                z = (X[n] + 1j * Y[n]) - (X[n + 1] + 1j * Y[n + 1])
                angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                         np.angle(z * np.exp(-1j * np.pi / 5.0))]
                X_arr = [X[n], X[n]]
                Y_arr = [Y[n], Y[n]]
                fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                             color=colours[id_colour], line_width=1, alpha=0.5)

                n = len(model_time_serie[i]['x']) - 2
                z = (X[n] + 1j * Y[n]) - (X[n + 1] + 1j * Y[n + 1])
                angle = [np.angle(z * np.exp(1j * np.pi / 5)),
                         np.angle(z * np.exp(-1j * np.pi / 5))]
                X_arr = [X[n + 1], X[n + 1]]
                Y_arr = [Y[n + 1], Y[n + 1]]
                fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                             color=colours[id_colour], line_width=1, alpha=0.5)

                XX = np.array(
                    X[rang_c] + 1j * Y[rang_c] + best_model['rho'] * np.exp(
                        1j * np.linspace(0, 2.0 * np.pi, 100)))
                # fig_curr.line(XX.real, XX.imag, line_width=0.5, color='black', alpha=0.5)
                fig_curr.patch(XX.real, XX.imag, color='black', alpha=0.3)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Trajectories (data)
            # ^^^^^^^^^^^^^^^^^^^
            fig_curr.circle('x', 'y', size=8, color='colour', alpha=0.5,
                            source=source)

            # Rotation for Earth + Spitzer
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if len(model_time_serie) == 2:
                # Find the vector D in (x, y) at t0
                source_t0_earth = np.array([ \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['x'],
                             kind='linear')(best_model['t0']), \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['y'],
                             kind='linear')(best_model['t0'])
                ])

                source_t0_earth_pente = np.array([ \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['x'],
                             kind='linear')(best_model['t0'] + 0.1), \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['y'],
                             kind='linear')(best_model['t0'] + 0.1)
                ])
                source_t0_earth_pente = (
                                            source_t0_earth_pente - source_t0_earth) / 0.1

                source_t0_spitzer = np.array([interp1d(
                    model_time_serie[1]['dates'],
                    model_time_serie[1]['x'], kind='linear')(best_model['t0']), \
                                              interp1d(
                                                  model_time_serie[1]['dates'],
                                                  model_time_serie[1]['y'],
                                                  kind='linear')(
                                                  best_model['t0'])])

                source_t0_spitzer_pente = np.array([ \
                    interp1d(model_time_serie[1]['dates'],
                             model_time_serie[1]['x'], kind='linear')(
                        best_model['t0'] + 0.1), \
                    interp1d(model_time_serie[1]['dates'],
                             model_time_serie[1]['y'], kind='linear')(
                        best_model['t0'] + 0.1)
                ])
                source_t0_spitzer_pente = (
                                              source_t0_spitzer_pente - source_t0_spitzer) / 0.1

                D_t0_xy = source_t0_spitzer - source_t0_earth

                # Angle between D in (x,y) and (E,N)
                # Caution: we rotate (x,y) by pi/2, so y is towards left, x towards
                # top. Now we can compute the rotation angle beetween \Delta\zeta
                # and D in (E,N). This angle + pi/2 gives the rotation angle to draw
                # trajectories in (E,N). All this is necessary because (E,N,T) is
                # equivalent to (y, x, T) with T is the target.
                D_xy_c = (D_t0_xy[0] + 1j * D_t0_xy[1]) * np.exp(1j * np.pi / 2.0)
                D_c = D_enm[0] + 1j * D_enm[1]

                alpha1 = np.angle(D_xy_c, deg=False)
                alpha2 = np.angle(D_c, deg=False)
                # epsilon = (angle_between(D_t0_xy, np.array([D_enm[0], D_enm[1]])))
                epsilon = (angle_between(np.array([D_xy_c.real, D_xy_c.imag]),
                                         np.array(
                                             [D_enm[0], D_enm[1]]))) + np.pi / 2.0
                rotation = np.exp(1j * epsilon)
                # print alpha1*180.0/np.pi, alpha2*180.0/np.pi, epsilon*180.0/np.pi

                # Unit vectors in xy
                x_hat_xy = 1.0
                y_hat_xy = 1j
                e_hat_xy = x_hat_xy * np.exp(1j * epsilon)
                n_hat_xy = y_hat_xy * np.exp(1j * epsilon)

                # Unit vectors in EN
                e_hat_en = 1.0
                n_hat_en = 1j
                x_hat_en = e_hat_en * np.exp(-1j * epsilon)
                y_hat_en = n_hat_en * np.exp(-1j * epsilon)

                # D in (x,y)
                palette = plt.get_cmap('Paired')
                id_palette = 0.090909  # from 0 to 11

                X = np.linspace(source_t0_earth[0], source_t0_spitzer[0], 2)
                Y = np.linspace(source_t0_earth[1], source_t0_spitzer[1], 2)
                n = 9

                fig_curr.line(X, Y, line_width=2, line_dash='dashed',
                              color='green', alpha=1)
                fig_curr.circle(X[0], Y[0], size=15, color='green', alpha=0.5)
                fig_curr.circle(X[1], Y[1], size=15, color='green', alpha=0.5)

                # ax_curr.plot(X, Y, dashes=(4, 2), lw=1,
                #                     color=palette(n * id_palette), alpha=1, zorder=20)
                # ax_curr.scatter(X[0], Y[0], 15, marker="o", linewidths=0.3,
                #                         facecolors=palette(n * id_palette), edgecolors='k',
                #                         alpha=0.8, zorder=20)
                # ax_curr.scatter(X[1], Y[1], 15, marker="o", linewidths=0.3,
                #                         facecolors=palette(n * id_palette), edgecolors='k',
                #                         alpha=0.8, zorder=20)
                #

            # Layout
            # ^^^^^^
            fig_curr.xaxis.axis_label = u'\u03B8\u2081 (Einstein units)'
            fig_curr.yaxis.axis_label = u'\u03B8\u2082 (Einstein units)'
            fig_curr.xaxis.axis_label_text_font = 'helvetica'
            fig_curr.yaxis.axis_label_text_font = 'helvetica'
            fig_curr.xaxis.axis_label_text_font_size = '10pt'
            fig_curr.yaxis.axis_label_text_font_size = '10pt'

            fig_curr.min_border_top = 10
            fig_curr.min_border_bottom = 0
            fig_curr.min_border_left = 0

            fig_curr.xgrid.grid_line_color = None
            fig_curr.ygrid.grid_line_color = None

            # ..................................................................
            #    Plot caustic 2 : pc2
            # ..................................................................

            # if 0:
            if len(model_time_serie) == 2:

                # Caution: we draw in (East, North) with East towards the left hand side.
                # So, X must be Eastern component, Y must be Northern component.
                # Then, always plot -X, Y to get plots in (West, North) frame.
                # Finaly reverse x-axis to get back to the East towards left.

                # Preparation
                # ^^^^^^^^^^^
                time_serie.update({'-x_complex': -(
                (time_serie['x'] + 1j * time_serie['y']) * rotation).real})
                time_serie.update({'y_complex': (
                (time_serie['x'] + 1j * time_serie['y']) * rotation).imag})

                # Plot
                # ^^^^
                source = ColumnDataSource(time_serie)

                hover_pc2 = HoverTool(
                    tooltips=[
                        ("ID", "@id{int}"),
                        ("Obs", "@obs"),
                        ("Date", "@dates{1.11}")
                    ]
                )

                xmin = -1.0
                xmax = 1.0
                ymin = -1.0
                ymax = 1.0
                tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap",
                         hover_pc2]
                fig = np.append(fig, \
                                bplt.figure(toolbar_location="above",
                                            plot_width=560, plot_height=600,
                                            x_range=(xmax, xmin),
                                            y_range=(ymin, ymax), \
                                            title=None, min_border=10,
                                            min_border_left=50, tools=tools))

                fig_curr = fig[3]

                # Caustic
                # ^^^^^^^
                if q > 1e-9:
                    caustic = caustic * rotation
                    fig_curr.circle(-caustic.real, caustic.imag, size=0.5, color='red',
                                    alpha=0.5)

                    X = (s0 * q / (1 + q)) * rotation
                    fig_curr.circle(-X.real, X.imag, size=5, color='orange', alpha=1)

                    X = (s0 * q / (1 + q) - s0) * rotation
                    fig_curr.circle(-X.real, X.imag, size=5, color='orange', alpha=1)
                else:
                    fig_curr.circle(0, 0, size=5, color='orange', alpha=1)

                # Trajectories
                # ^^^^^^^^^^^^
                colours = ['black', '#297CC4']
                id_colour = 0
                for i in xrange(len(locations)):
                    X = (model_time_serie[i]['x'] + 1j * model_time_serie[i][
                        'y']) * rotation
                    fig_curr.line(-X.real, X.imag, line_width=2,
                                  color=colours[id_colour], alpha=1)

                    # Arrows
                    n = 0
                    z = X[n] - X[n + 1]
                    angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                             np.angle(z * np.exp(-1j * np.pi / 5.0))]
                    X_arr = [-X[n].real, -X[n].real]
                    Y_arr = [X[n].imag, X[n].imag]
                    fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                                 color=colours[id_colour], line_width=1, alpha=0.5)

                    n = len(model_time_serie[i]['x']) - 2
                    z = X[n] - X[n + 1]
                    angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                             np.angle(z * np.exp(-1j * np.pi / 5.0))]
                    X_arr = [-X[n + 1].real, -X[n + 1].real]
                    Y_arr = [X[n + 1].imag, X[n + 1].imag]
                    fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                                 color=colours[id_colour], line_width=1, alpha=0.5)

                    if id_colour < len(colours) - 1:
                        id_colour = id_colour + 1
                    else:
                        id_colour = 0

                # Trajectories (data)
                # ^^^^^^^^^^^^^^^^^^^
                fig_curr.circle('-x_complex', 'y_complex', size=8, color='colour',
                                alpha=0.5, source=source)

                # Some specific positions
                # ^^^^^^^^^^^^^^^^^^^^^^^
                # D in (e,n) at tO
                A = (source_t0_earth[0] + 1j * source_t0_earth[1]) * rotation
                B = (source_t0_spitzer[0] + 1j * source_t0_spitzer[1]) * rotation
                X = np.linspace(A.real, B.real, 2)
                Y = np.linspace(A.imag, B.imag, 2)
                n = 9
                fig_curr.line(-X, Y, line_width=2, line_dash='dashed',
                              color='green', alpha=1)
                fig_curr.circle(-X[0], Y[0], size=15, color='green', alpha=0.5)
                fig_curr.circle(-X[1], Y[1], size=15, color='green', alpha=0.5)

                # Layout
                # ^^^^^^
                fig_curr.xaxis.axis_label = 'theta / theta_E (East)'
                fig_curr.yaxis.axis_label = 'theta / theta_E (North)'
                fig_curr.xaxis.axis_label_text_font = 'helvetica'
                fig_curr.yaxis.axis_label_text_font = 'helvetica'
                fig_curr.xaxis.axis_label_text_font_size = '10pt'
                fig_curr.yaxis.axis_label_text_font_size = '10pt'

                fig_curr.min_border_top = 10
                fig_curr.min_border_bottom = 0
                fig_curr.min_border_left = 0

                fig_curr.xgrid.grid_line_color = None
                fig_curr.ygrid.grid_line_color = None

            # ------------------------------------------------------------------
            #   Save the html page
            # ------------------------------------------------------------------
            if len(model_time_serie) == 2:
                final = blyt.column(fig[0], fig[1], blyt.row(fig[2], fig[3]))
                bplt.save(final)
            if len(model_time_serie) != 2:
                final = blyt.column(fig[0], fig[1], fig[2])
                # final = blyt.column(fig[0], fig[1])
                bplt.save(final)

            # ------------------------------------------------------------------
            #   Modify the html page
            # ------------------------------------------------------------------
            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                'RelativePaths', 'Plots')

            if cfgsetup.getboolean('Plotting', 'Data'):
                filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                    'RelativePaths', 'Plots')

                if (best_model['fullid'] == -1) & flag_fix:
                    filename = filename + cfgsetup.get('Controls',
                                                       'Archive') + '-summary.html'
                    title = cfgsetup.get('EventDescription',
                                         'Name') + ': best model last MCMC'
                elif flag_fix:
                    filename = filename + cfgsetup.get('Controls',
                                                       'Archive') + '-summary-' \
                               + repr(best_model['fullid']) + '.html'
                    title = cfgsetup.get('EventDescription',
                                         'Name') + ' - Model # ' + repr(
                        best_model['fullid'])
                else:
                    filename = filename + cfgsetup.get('Controls',
                                                       'Archive') + '-summary-fix.html'
                    title = cfgsetup.get('EventDescription',
                                         'Name') + ' - Model fix'






            file = open(filename, 'r')
            file_new = ''
            for line in file:
                # print line.strip()[:7]
                if line.strip()[:7] == '<title>':
                    file_new = file_new \
                            + '        <style type="text/css">\n' \
                            + '        p.p1 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 43.0px; font: 36.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n'\
                            + '        p.p2 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 21.0px; font: 18.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n'\
                            + '        p.p3 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 15.0px; font: 12.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n'\
                            + '        p.p4 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000; min-height: 17.0px}\n'\
                            + '        p.p5 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n'\
                            + '        p.p6 {margin: 0.0px 0.0px 12.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n'\
                            + '        p.p7 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n'\
                            + '        span.s1 {font-kerning: none}\n'\
                            + '        span.s10 {font: 14.0px "Lucida Grande"; color: #585858}\n'\
                            + '        hr {\n'\
                            + '        display: block;\n'\
                            + '        margin-top: 0.5em;\n'\
                            + '        margin-bottom: 0.5em;\n'\
                            + '        margin-left: auto;\n'\
                            + '        margin-right: auto;\n'\
                            + '        border-style: inset;\n'\
                            + '        border-width: 1px;\n'\
                            + '        }\n'\
                            + '        </style>\n'\
                            + '        <title>' + 'muLAn ' + cfgsetup.get('EventDescription', 'Name')[4:] + '/' + cfgsetup.get('Controls', 'Archive') + '#' + repr(best_model['fullid']) + '</title>\n'\
                            + '        <meta name="Author" content="Clement Ranc">\n'
                elif line.strip()[:7] == '</head>':
                    file_new = file_new\
                               + '        <script type="text/x-mathjax-config">\n'\
                               + "        MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\\(','\\\)']]}});\n"\
                               + '        </script>\n'\
                               + '        <script type="text/javascript" async  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML"></script>\n'\
                               + '    </head>\n'
                elif line.strip()[:6] == '<body>':
                    file_new = file_new \
                               + '    <body>\n\n' \
                               + '<p class="p1"><span class="s1"><b>' + title + '</b></span></p>\n' \
                               + '<p class="p2"><span class="s1"><br>\n' \
                               + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(t_0\)    = ' + repr(best_model['t0']) + ' days</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(u_0\)       = ' + repr(best_model['u0']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(t_\mathrm{E}\)       = ' + repr(best_model['tE']) + ' days</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\\rho\)      = ' + repr(best_model['rho']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\pi_\mathrm{EN}\)     = ' + repr(best_model['piEN']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\pi_\mathrm{EE}\)     = ' + repr(best_model['piEE']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(s\)        = ' + repr(best_model['s']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(q\)        = ' + repr(best_model['q']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\\alpha\)    = ' + repr(best_model['alpha']) + ' radians</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\mathrm{d}\\alpha/\mathrm{d}t\)= ' + repr(best_model['dalpha']) + ' radians/years</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\mathrm{d}s/\mathrm{d}t\)    = ' + repr(best_model['ds']) + ' years^-1</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\chi^2\)     = ' + repr(chi2_flux) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">\(\chi^2/\mathrm{dof}\) = ' + repr(chi2dof_flux) + '</span></p>\n' \
                               + '<p class="p2"><span class="s1"><br>\n' \
                               + '</span></p>\n'
                elif line.strip()[:7] == '</body>':
                    file_new = file_new \
                               + '        <BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n' \
                               + '        <BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n' \
                               + '        <hr>\n' \
                               + '        <BR>\n' \
                               + '        <footer>\n'\
                               + '        <p class="p7"><span class="s10">Modelling and page by muLAn (MicroLensing Analysis software).</span></p>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        </footer>\n' \
                               + '        </body>\n'
                else:
                    file_new = file_new + line
            file.close()

            file = open(filename, 'w')
            file.write(file_new)
            file.close()



            if 0:
                # PLOT STATISTIC RESIDUAL IN MAG ==================================

                for j in xrange(len(observatories)):
                    cond2 = (time_serie['obs'] == observatories[j])

                    # Configuration
                    # ------------------------------------------------------------------
                    moteur = 'defaultSmall_pdf'
                    # ------------------------------------------------------------------
                    # Conversions
                    in2cm = 2.54
                    cm2in = 1.0 / in2cm
                    in2pt = 72.27
                    pt2in = 1.0 / 72.27
                    cm2pt = cm2in * in2pt
                    pt2cm = 1.0 / cm2pt
                    # ------------------------------------------------------------------
                    # Configuration à compléter
                    fig_width_cm = 6.5
                    fig_height_cm = 4.5
                    path = cfgsetup.get('FullPaths', 'Event')\
                           + cfgsetup.get('RelativePaths', 'Plots')
                    filename = path + cfgsetup.get('Controls', 'Archive')\
                               + "-" + observatories[j]\
                               + '-Residuals_Statistics'
                    # ------------------------------------------------------------------
                    # Calculs pour configuration
                    fig_width_pt = fig_width_cm * cm2pt
                    fig_height_pt = fig_height_cm * cm2pt  # golden_mean = (math.sqrt(5)-1.0)/2.0 = 0.62
                    extension = moteur[-3:]
                    filename_moteur = full_path + 'plotconfig/matplotlibrc_' + moteur
                    pylab.rcParams.update(
                        mpl.rc_params_from_file(filename_moteur, fail_on_error=True))
                    fig_size = np.array([fig_width_pt, fig_height_pt]) * pt2in
                    # ------------------------------------------------------------------

                    fig1 = plt.figure('Figure', figsize=fig_size)

                    # ..................................................................
                    #   Plot 1
                    # ..................................................................
                    layout = [1.0, 0.8, 5, 3.1]  # en cm
                    kde = np.array(layout) * cm2in / np.array(
                        [fig_size[0], fig_size[1], fig_size[0], fig_size[1]])
                    ax_curr = fig1.add_axes(kde)

                    grandeur = time_serie['residus'][cond2]
                    nb_bins = 2*int(np.sqrt(np.max([len(grandeur), 3])))

                    lim_stat = [np.min(grandeur), np.max(grandeur)]

                    X = grandeur
                    hist, bin_edges = np.histogram(X, bins=nb_bins,
                                                   range=(lim_stat), density=1)
                    hist_plot = np.append(hist, [hist[-1]]) / np.max(hist)
                    ax_curr.step(bin_edges, hist_plot, 'k-', where="mid",
                                 zorder=2)

                    # Model
                    model = GaussianMixture(n_components=1).fit(np.atleast_2d(grandeur).T)
                    X = np.atleast_2d(np.linspace(lim_stat[0], lim_stat[1], 1000)).T
                    logprob = model.score_samples(X)
                    pdf_individual = np.exp(logprob)
                    pdf_individual = pdf_individual / np.max(pdf_individual)

                    ax_curr.plot(X.T[0], pdf_individual, ls='-', color='b', zorder=1)

                    mean = np.mean(X.T[0])
                    rms = np.std(X.T[0])

                    chat = "mean {:.4f}\nstd dev {:.4f}".format(mean, rms)
                    position = [1, 1]
                    ax_curr.annotate(chat, xy=position, xycoords='axes fraction',
                                     ha="right", va="center", color='k', fontsize=5,
                                     backgroundcolor='w', zorder=100)

                    # Limits
                    # ^^^^^^
                    # ax_curr.set_xlim(0, 3.5)
                    ax_curr.set_ylim(0, 1.05)
                    # ax_curr.set_ylim(ax_curr.get_ylim()[1], ax_curr.get_ylim()[0])

                    # ax_curr.set_xscale('log')
                    # ax_curr.set_yscale('log')

                    # Ticks
                    # ^^^^^
                    # ax_curr.xaxis.set_major_locator(MultipleLocator(0.4))
                    ax_curr.xaxis.set_major_locator(MaxNLocator(5))
                    minor = 0.5 * (np.roll(ax_curr.get_xticks(), -1) - ax_curr.get_xticks())[0]
                    minor_locator = MultipleLocator(minor)
                    ax_curr.xaxis.set_minor_locator(minor_locator)

                    ax_curr.yaxis.set_major_locator(MultipleLocator(0.2))
                    # ax_curr.yaxis.set_major_locator(MaxNLocator(4))
                    minor = 0.5 * (np.roll(ax_curr.get_yticks(), -1) - ax_curr.get_yticks())[0]
                    minorLocator = MultipleLocator(minor)
                    ax_curr.yaxis.set_minor_locator(minorLocator)

                    # Legend
                    # ^^^^^^
                    ax_curr.set_xlabel(ur"%s" % ("$\sigma$ (mag)"), labelpad=0)
                    ax_curr.set_ylabel(ur"%s" % ("count"), labelpad=3)

                    # # --> Options de fin
                    #
                    # for tick in ax_curr.get_yaxis().get_major_ticks():
                    #     tick.set_pad(2)
                    #     tick.label1 = tick._get_text1()
                    # for tick in ax_curr.get_xaxis().get_major_ticks():
                    #     tick.set_pad(2)
                    #     tick.label1 = tick._get_text1()

                    # ..................................................................
                    #   SAVE FIGURE
                    # ..................................................................
                    if 0:
                        plt.show()
                    fig1.savefig(filename + "." + extension, transparent=False,
                                 dpi=1400)
                    plt.close()





                # =================================================================












        else:
            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                'RelativePaths', 'Plots')

            if (best_model['fullid'] == -1) & flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ': best model last MCMC'
            elif flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification-' \
                           + repr(best_model['fullid']) + '.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model # ' + repr(
                    best_model['fullid'])
            else:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification-fix.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model fix'

            if os.path.exists(filename):
                os.remove(filename)
            bplt.output_file(filename)
            fig = np.array([])
            plot_counter = 0

            observatories = np.unique(time_serie['obs'])

            # ..................................................................
            #    Plot light curve : amplification
            # ..................................................................
            tmin = float(options.split('/')[0].split('-')[0].strip())
            tmax = float(options.split('/')[0].split('-')[1].strip())

            ymin = 0.9
            ymax = 0
            for i in xrange(len(locations)):
                ymax_temp = np.max(model_time_serie[i]['amp'])
                ymax = np.max(np.array([ymax, ymax_temp]))

            tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap"]
            fig = np.append(fig, \
                            bplt.figure(toolbar_location="above", plot_width=1200,
                                        plot_height=600, x_range=(tmin, tmax),
                                        y_range=(ymin, ymax), \
                                        title=None, min_border=10,
                                        min_border_left=50, tools=tools))

            fig_curr = fig[plot_counter]

            # Annotations
            # ^^^^^^^^^^^
            colours = ['black', '#297CC4']
            id_colour = 0
            for i in xrange(len(locations)):
                name = 'Models_' + locations[i]
                models_temp = model_param[name]
                name = 'DateRanges_' + locations[i]
                dates_temp = model_param[name]

                for j in xrange(len(models_temp)):
                    tmin = float((dates_temp[j]).split('-')[0].strip())
                    tmax = float((dates_temp[j]).split('-')[1].strip())

                    X = np.ones(2) * tmin
                    Y = np.linspace(-100000, 100000, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                    X = np.ones(2) * tmax
                    Y = np.linspace(-100000, 100000, 2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                    X = np.linspace(-100000, 100000, 2)
                    Y = 1.0 * np.ones(2)
                    fig_curr.line(X, Y, line_width=0.5, line_dash='dashed',
                                  color=colours[id_colour], alpha=0.4)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Amplification models
            # ^^^^^^^^^^^^^^^^^^^^
            colours = ['black', '#297CC4', 'green']
            id_colour = 0
            for i in xrange(len(locations)):
                X = model_time_serie[i]['dates']
                Y = model_time_serie[i]['amp']
                fig_curr.line(X, Y, line_width=2, color=colours[id_colour],
                              alpha=1)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Layout
            # ^^^^^^
            fig_curr.xaxis.axis_label = 'HJD - 2,450,000'
            fig_curr.yaxis.axis_label = 'Amplification'
            fig_curr.xaxis.axis_label_text_font = 'helvetica'
            fig_curr.yaxis.axis_label_text_font = 'helvetica'
            fig_curr.xaxis.axis_label_text_font_size = '10pt'
            fig_curr.yaxis.axis_label_text_font_size = '10pt'

            fig_curr.min_border_top = 10
            fig_curr.min_border_bottom = 0
            fig_curr.min_border_left = 0

            fig_curr.xgrid.grid_line_color = None
            fig_curr.ygrid.grid_line_color = None

            plot_counter = plot_counter + 1

            # ..................................................................
            #    Plot caustic 1 : pc1
            # ..................................................................

            # Plot
            # ^^^^
            xmin = -1.0
            xmax = 1.0
            ymin = -1.0
            ymax = 1.0
            tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap"]
            fig = np.append(fig, \
                            bplt.figure(toolbar_location="above", plot_width=600,
                                        plot_height=560, x_range=(xmin, xmax),
                                        y_range=(ymin, ymax), \
                                        title=None, min_border=10,
                                        min_border_left=50, tools=tools))

            fig_curr = fig[plot_counter]

            # Caustic
            # ^^^^^^^
            # Case of lens orbital rotation
            try:
                time_caustic = options.split('/')[2].replace('[','').replace(']','').split('-')
                time_caustic = np.array([float(a.strip()) for a in time_caustic])
                n_caustics = len(time_caustic)
                nb_pts_caus = 1000
                alpha, s = lens_rotation(alpha0, s0, dalpha, ds, time_caustic, tb)
                color_caustics = np.array(['Orange', 'SeaGreen', 'LightSeaGreen', 'CornflowerBlue', 'DarkViolet'])
            except:
                nb_pts_caus = 1000
                n_caustics = 0

            if q > 1e-9:
                if n_caustics > 0:
                    for i in xrange(n_caustics):
                        # > Courbes critiques et caustiques
                        delta = 2.0 * np.pi / (nb_pts_caus - 1.0)
                        phi = np.arange(nb_pts_caus) * delta
                        critic = np.array([])
                        for phi_temp in phi[:]:
                            critic = np.append(critic, critic_roots(s[i], q, phi_temp))
                        caustic = critic - 1 / (1 + q) * (1 / critic.conjugate() + q / (critic.conjugate() + s[i]))
                        caustic = caustic + GL1  # From Cassan (2008) to CM
                        caustic = caustic * np.exp(1j*(alpha[i]-alpha0))

                        fig_curr.circle(caustic.real, caustic.imag, size=0.5, color=color_caustics[0], alpha=0.5)
                        color_caustics = np.roll(color_caustics, -1)

                # > Courbes critiques et caustiques
                delta = 2.0 * np.pi / (nb_pts_caus - 1.0)
                phi = np.arange(nb_pts_caus) * delta
                critic = np.array([])
                for phi_temp in phi[:]:
                    critic = np.append(critic, critic_roots(s0, q, phi_temp))
                caustic = critic - 1 / (1 + q) * (
                1 / critic.conjugate() + q / (critic.conjugate() + s0))
                caustic = caustic + GL1  # From Cassan (2008) to CM

                fig_curr.circle(caustic.real, caustic.imag, size=0.5, color='red',
                                alpha=0.5)
                fig_curr.circle(GL1.real, GL1.imag, size=5, color='orange', alpha=1)
                fig_curr.circle(GL2.real, GL2.imag, size=5, color='orange', alpha=1)
            else:
                fig_curr.circle(0, 0, size=5, color='orange', alpha=1)

            # Trajectories
            # ^^^^^^^^^^^^
            colours = ['black', '#297CC4']
            id_colour = 0

            for i in xrange(len(locations)):
                temp = np.array([abs(a - best_model['t0']) for a in
                                 model_time_serie[i]['dates']])
                rang_c = np.where(temp == np.min(temp))[0][0]

                X = model_time_serie[i]['x']
                Y = model_time_serie[i]['y']
                fig_curr.line(X, Y, line_width=2, color=colours[id_colour],
                              alpha=1)

                # Arrows
                n = 0
                z = (X[n] + 1j * Y[n]) - (X[n + 1] + 1j * Y[n + 1])
                angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                         np.angle(z * np.exp(-1j * np.pi / 5.0))]
                X_arr = [X[n], X[n]]
                Y_arr = [Y[n], Y[n]]
                fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                             color=colours[id_colour], line_width=1, alpha=0.5)

                n = len(model_time_serie[i]['x']) - 2
                z = (X[n] + 1j * Y[n]) - (X[n + 1] + 1j * Y[n + 1])
                angle = [np.angle(z * np.exp(1j * np.pi / 5)),
                         np.angle(z * np.exp(-1j * np.pi / 5))]
                X_arr = [X[n + 1], X[n + 1]]
                Y_arr = [Y[n + 1], Y[n + 1]]
                fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                             color=colours[id_colour], line_width=1, alpha=0.5)

                XX = np.array(
                    X[rang_c] + 1j * Y[rang_c] + best_model['rho'] * np.exp(
                        1j * np.linspace(0, 2.0 * np.pi, 100)))
                # fig_curr.line(XX.real, XX.imag, line_width=0.5, color='black', alpha=0.5)
                fig_curr.patch(XX.real, XX.imag, color='black', alpha=0.3)

                if id_colour < len(colours) - 1:
                    id_colour = id_colour + 1
                else:
                    id_colour = 0

            # Rotation for Earth + Spitzer
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if len(model_time_serie) == 2:
                # Find the vector D in (x, y) at t0
                source_t0_earth = np.array([ \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['x'],
                             kind='linear')(best_model['t0']), \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['y'],
                             kind='linear')(best_model['t0'])
                ])

                source_t0_earth_pente = np.array([ \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['x'],
                             kind='linear')(best_model['t0'] + 0.1), \
                    interp1d(model_time_serie[0]['dates'],
                             model_time_serie[0]['y'],
                             kind='linear')(best_model['t0'] + 0.1)
                ])
                source_t0_earth_pente = (
                                            source_t0_earth_pente - source_t0_earth) / 0.1

                source_t0_spitzer = np.array([interp1d(
                    model_time_serie[1]['dates'],
                    model_time_serie[1]['x'], kind='linear')(best_model['t0']), \
                                              interp1d(
                                                  model_time_serie[1]['dates'],
                                                  model_time_serie[1]['y'],
                                                  kind='linear')(
                                                  best_model['t0'])])

                source_t0_spitzer_pente = np.array([ \
                    interp1d(model_time_serie[1]['dates'],
                             model_time_serie[1]['x'], kind='linear')(
                        best_model['t0'] + 0.1), \
                    interp1d(model_time_serie[1]['dates'],
                             model_time_serie[1]['y'], kind='linear')(
                        best_model['t0'] + 0.1)
                ])
                source_t0_spitzer_pente = (
                                              source_t0_spitzer_pente - source_t0_spitzer) / 0.1

                D_t0_xy = source_t0_spitzer - source_t0_earth

                # Angle between D in (x,y) and (E,N)
                # Caution: we rotate (x,y) by pi/2, so y is towards left, x towards
                # top. Now we can compute the rotation angle beetween \Delta\zeta
                # and D in (E,N). This angle + pi/2 gives the rotation angle to draw
                # trajectories in (E,N). All this is necessary because (E,N,T) is
                # equivalent to (y, x, T) with T is the target.
                D_xy_c = (D_t0_xy[0] + 1j * D_t0_xy[1]) * np.exp(1j * np.pi / 2.0)
                D_c = D_enm[0] + 1j * D_enm[1]

                alpha1 = np.angle(D_xy_c, deg=False)
                alpha2 = np.angle(D_c, deg=False)
                # epsilon = (angle_between(D_t0_xy, np.array([D_enm[0], D_enm[1]])))
                epsilon = (angle_between(np.array([D_xy_c.real, D_xy_c.imag]),
                                         np.array(
                                             [D_enm[0], D_enm[1]]))) + np.pi / 2.0
                rotation = np.exp(1j * epsilon)
                # print alpha1*180.0/np.pi, alpha2*180.0/np.pi, epsilon*180.0/np.pi

                # Unit vectors in xy
                x_hat_xy = 1.0
                y_hat_xy = 1j
                e_hat_xy = x_hat_xy * np.exp(1j * epsilon)
                n_hat_xy = y_hat_xy * np.exp(1j * epsilon)

                # Unit vectors in EN
                e_hat_en = 1.0
                n_hat_en = 1j
                x_hat_en = e_hat_en * np.exp(-1j * epsilon)
                y_hat_en = n_hat_en * np.exp(-1j * epsilon)

                # D in (x,y)
                palette = plt.get_cmap('Paired')
                id_palette = 0.090909  # from 0 to 11

                X = np.linspace(source_t0_earth[0], source_t0_spitzer[0], 2)
                Y = np.linspace(source_t0_earth[1], source_t0_spitzer[1], 2)
                n = 9

                fig_curr.line(X, Y, line_width=2, line_dash='dashed',
                              color='green', alpha=1)
                fig_curr.circle(X[0], Y[0], size=15, color='green', alpha=0.5)
                fig_curr.circle(X[1], Y[1], size=15, color='green', alpha=0.5)

                # ax_curr.plot(X, Y, dashes=(4, 2), lw=1,
                #                     color=palette(n * id_palette), alpha=1, zorder=20)
                # ax_curr.scatter(X[0], Y[0], 15, marker="o", linewidths=0.3,
                #                         facecolors=palette(n * id_palette), edgecolors='k',
                #                         alpha=0.8, zorder=20)
                # ax_curr.scatter(X[1], Y[1], 15, marker="o", linewidths=0.3,
                #                         facecolors=palette(n * id_palette), edgecolors='k',
                #                         alpha=0.8, zorder=20)
                #

            # Layout
            # ^^^^^^
            fig_curr.xaxis.axis_label = 'theta_x / theta_E'
            fig_curr.yaxis.axis_label = 'theta_y / theta_E'
            fig_curr.xaxis.axis_label_text_font = 'helvetica'
            fig_curr.yaxis.axis_label_text_font = 'helvetica'
            fig_curr.xaxis.axis_label_text_font_size = '10pt'
            fig_curr.yaxis.axis_label_text_font_size = '10pt'

            fig_curr.min_border_top = 10
            fig_curr.min_border_bottom = 0
            fig_curr.min_border_left = 0

            fig_curr.xgrid.grid_line_color = None
            fig_curr.ygrid.grid_line_color = None

            plot_counter = plot_counter + 1

            # ..................................................................
            #    Plot caustic 2 : pc2
            # ..................................................................

            if len(model_time_serie) == 2:

                # Caution: we draw in (East, North) with East towards the left hand side.
                # So, X must be Eastern component, Y must be Northern component.
                # Then, always plot -X, Y to get plots in (West, North) frame.
                # Finaly reverse x-axis to get back to the East towards left.

                # Plot
                # ^^^^
                xmin = -1.0
                xmax = 1.0
                ymin = -1.0
                ymax = 1.0
                tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap"]
                fig = np.append(fig, \
                                bplt.figure(toolbar_location="above",
                                            plot_width=600, plot_height=560,
                                            x_range=(xmax, xmin),
                                            y_range=(ymin, ymax), \
                                            title=None, min_border=10,
                                            min_border_left=50, tools=tools))

                fig_curr = fig[plot_counter]

                # Caustic
                # ^^^^^^^
                if q > 1e-9:
                    caustic = caustic * rotation
                    fig_curr.circle(-caustic.real, caustic.imag, size=0.5, color='red',
                                    alpha=0.5)

                    X = (s0 * q / (1 + q)) * rotation
                    fig_curr.circle(-X.real, X.imag, size=5, color='orange', alpha=1)

                    X = (s0 * q / (1 + q) - s) * rotation
                    fig_curr.circle(-X.real, X.imag, size=5, color='orange', alpha=1)
                else:
                    fig_curr.circle(0, 0, size=5, color='orange', alpha=1)

                # Trajectories
                # ^^^^^^^^^^^^
                colours = ['black', '#297CC4']
                id_colour = 0
                for i in xrange(len(locations)):
                    X = (model_time_serie[i]['x'] + 1j * model_time_serie[i][
                        'y']) * rotation
                    fig_curr.line(-X.real, X.imag, line_width=2,
                                  color=colours[id_colour], alpha=1)

                    # Arrows
                    n = 0
                    z = X[n] - X[n + 1]
                    angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                             np.angle(z * np.exp(-1j * np.pi / 5.0))]
                    X_arr = [-X[n].real, -X[n].real]
                    Y_arr = [X[n].imag, X[n].imag]
                    fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                                 color=colours[id_colour], line_width=1, alpha=0.5)

                    n = len(model_time_serie[i]['x']) - 2
                    z = X[n] - X[n + 1]
                    angle = [np.angle(z * np.exp(1j * np.pi / 5.0)),
                             np.angle(z * np.exp(-1j * np.pi / 5.0))]
                    X_arr = [-X[n + 1].real, -X[n + 1].real]
                    Y_arr = [X[n + 1].imag, X[n + 1].imag]
                    fig_curr.ray(X_arr, Y_arr, length=0.05, angle=angle,
                                 color=colours[id_colour], line_width=1, alpha=0.5)

                    if id_colour < len(colours) - 1:
                        id_colour = id_colour + 1
                    else:
                        id_colour = 0

                # Some specific positions
                # ^^^^^^^^^^^^^^^^^^^^^^^
                # D in (e,n) at tO
                A = (source_t0_earth[0] + 1j * source_t0_earth[1]) * rotation
                B = (source_t0_spitzer[0] + 1j * source_t0_spitzer[1]) * rotation
                X = np.linspace(A.real, B.real, 2)
                Y = np.linspace(A.imag, B.imag, 2)
                n = 9
                fig_curr.line(-X, Y, line_width=2, line_dash='dashed',
                              color='green', alpha=1)
                fig_curr.circle(-X[0], Y[0], size=15, color='green', alpha=0.5)
                fig_curr.circle(-X[1], Y[1], size=15, color='green', alpha=0.5)

                # Layout
                # ^^^^^^
                fig_curr.xaxis.axis_label = 'theta / theta_E (East)'
                fig_curr.yaxis.axis_label = 'theta / theta_E (North)'
                fig_curr.xaxis.axis_label_text_font = 'helvetica'
                fig_curr.yaxis.axis_label_text_font = 'helvetica'
                fig_curr.xaxis.axis_label_text_font_size = '10pt'
                fig_curr.yaxis.axis_label_text_font_size = '10pt'

                fig_curr.min_border_top = 10
                fig_curr.min_border_bottom = 0
                fig_curr.min_border_left = 0

                fig_curr.xgrid.grid_line_color = None
                fig_curr.ygrid.grid_line_color = None

                plot_counter = plot_counter + 1

            # ------------------------------------------------------------------
            #   Save the html page
            # ------------------------------------------------------------------
            if len(model_time_serie) == 2:
                final = blyt.column(fig[0], blyt.row(fig[1], fig[2]))
                bplt.save(final)
            if len(model_time_serie) != 2:
                final = blyt.column(fig[0], fig[1])
                # final = bplt.vplot(fig[2])
                bplt.save(final)

            # ------------------------------------------------------------------
            #   Modify the html page
            # ------------------------------------------------------------------
            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                'RelativePaths', 'Plots')

            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get(
                'RelativePaths', 'Plots')

            if (best_model['fullid'] == -1) & flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ': best model last MCMC'
            elif flag_fix:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification-' \
                           + repr(best_model['fullid']) + '.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model # ' + repr(
                    best_model['fullid'])
            else:
                filename = filename + cfgsetup.get('Controls',
                                                   'Archive') + '-amplification-fix.html'
                title = cfgsetup.get('EventDescription',
                                     'Name') + ' - Model fix'

            file = open(filename, 'r')
            file_new = ''
            for line in file:
                # print line.strip()[:7]
                if line.strip()[:7] == '<title>':
                    file_new = file_new \
                               + '        <style type="text/css">\n' \
                               + '        p.p1 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 43.0px; font: 36.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                               + '        p.p2 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 21.0px; font: 18.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                               + '        p.p3 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 15.0px; font: 12.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                               + '        p.p4 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000; min-height: 17.0px}\n' \
                               + '        p.p5 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n' \
                               + '        p.p6 {margin: 0.0px 0.0px 12.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n' \
                               + '        p.p7 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                               + '        span.s1 {font-kerning: none}\n' \
                               + '        span.s10 {font: 14.0px "Lucida Grande"; color: #585858}\n' \
                               + '        hr {\n' \
                               + '        display: block;\n' \
                               + '        margin-top: 0.5em;\n' \
                               + '        margin-bottom: 0.5em;\n' \
                               + '        margin-left: auto;\n' \
                               + '        margin-right: auto;\n' \
                               + '        border-style: inset;\n' \
                               + '        border-width: 1px;\n' \
                               + '        }\n' \
                               + '        </style>\n' \
                               + '        <title>' + 'muLAn ' + cfgsetup.get('EventDescription', 'Name')[4:] + '/' + cfgsetup.get('Controls', 'Archive') + '#'\
                                        + repr(best_model['fullid']) + '</title>\n' \
                               + '        <meta name="Author" content="Clement Ranc">\n'
                elif line.strip()[:6] == '<body>':
                    file_new = file_new \
                               + '    <body>\n\n' \
                               + '<p class="p1"><span class="s1"><b>' + title + '</b></span></p>\n' \
                               + '<p class="p2"><span class="s1"><br>\n' \
                               + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">t0       = ' + repr(
                        best_model['t0']) + ' days</span></p>\n' \
                               + '<p class="p3"><span class="s1">u0       = ' + repr(
                        best_model['u0']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">tE       = ' + repr(
                        best_model['tE']) + ' days</span></p>\n' \
                               + '<p class="p3"><span class="s1">rho      = ' + repr(
                        best_model['rho']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">piEN     = ' + repr(
                        best_model['piEN']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">piEE     = ' + repr(
                        best_model['piEE']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">s        = ' + repr(
                        best_model['s']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">q        = ' + repr(
                        best_model['q']) + '</span></p>\n' \
                               + '<p class="p3"><span class="s1">alpha    = ' + repr(
                        best_model['alpha']) + ' radians</span></p>\n' \
                               + '<p class="p3"><span class="s1">dalpha/dt= ' + repr(
                        best_model['dalpha']) + ' radians/years</span></p>\n'\
                               + '<p class="p3"><span class="s1">ds/dt    = ' + repr(
                        best_model['ds']) + ' years^-1</span></p>\n' \
                               + '<p class="p3"><span class="s1">chi2     = ' + repr(
                        best_model['chi2']) + ' radians</span></p>\n' \
                               + '<p class="p3"><span class="s1">chi2/dof = ' + repr(
                        best_model['chi2/dof']) + ' radians</span></p>\n' \
                               + '<p class="p2"><span class="s1"><br>\n' \
                               + '</span></p>\n'
                elif line.strip()[:7] == '</body>':
                    file_new = file_new \
                               + '        <BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n<BR>\n' \
                               + '        <hr>\n' \
                               + '        <BR>\n' \
                               + '        <footer>\n' \
                               + '        <p class="p7"><span class="s10">Modelling and page by muLAn (MicroLensing Analysis software).</span></p>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        <BR>\n' \
                               + '        </footer>\n' \
                               + '        </body>\n'
                else:
                    file_new = file_new + line
            file.close()

            file = open(filename, 'w')
            file.write(file_new)
            file.close()

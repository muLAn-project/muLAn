# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------

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
while i < len(text)-1:
   a = a + text[i] + '/'
   i = i + 1
full_path = a

filename = full_path + '../' + '.pythonexternallibpath'
file = open(filename, 'r')
for line in file:
    path_lib_ext=line
file.close()
if path_lib_ext != 'None':
    sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#   Packages
# ----------------------------------------------------------------------
import numpy as np
import pickle
from scipy.interpolate import interp1d
import pandas as pd
pd.set_option('display.width', 50)

# ==============================================================================
#   MAIN
# ==============================================================================
#def magnifcalc(timeserie, physical_properties, Ds=0, tb=None):
def magnifcalc(t, param, Ds=None, tb=None):

    t0 = param['t0']
    u0 = param['u0']
    tE = param['tE']
    rho = param['rho']
    gamma = param['gamma']
    piEN = param['piEN']
    piEE = param['piEE']

    if tb == None: tb = t0
    if Ds == None:
        Ds = dict({'N' : np.zeros(len(t)),\
                   'E' : np.zeros(len(t))})

    # Correction of the trajectory due to parallax
    DsN = Ds['N']
    DsE = Ds['E']

    fn_bo = full_path + "espl_inputs/tab_Bo.p"
    fn_b1 = full_path + "espl_inputs/tab_B1.p"

    return amplification_limb(t0, tE, u0, rho, gamma, t, piEN, piEE, DsN, DsE,
           fn_bo, fn_b1)

def f_bo_extrapol(x, x_min, x_max, x_tau, x_beta, x_u, f_bo, f_b1):

    """x should be an array."""
    cond_inside = np.logical_and(x <= x_max, x >= x_min)
    cond_above = np.where(x > x_max)
    cond_below = np.where(x < x_min)
    nb_inside = np.sum(cond_inside)
    nb_above = np.sum(x > x_max)
    nb_below = np.sum(x < x_min)
    data = pd.DataFrame({'z': np.array([]),
                         'bo': np.array([]),
                         'b1': np.array([]),
                         'tau': np.array([]),
                         'beta': np.array([]),
                         'u': np.array([]),
                         'es': np.array([])})
    if nb_inside > 0:
        data_inside = pd.DataFrame({'z': x[cond_inside],
                                    'tau': x_tau[cond_inside],
                                    'beta': x_beta[cond_inside],
                                    'u': x_u[cond_inside],
                                    'b1': f_b1(x[cond_inside]),
                                    'bo': f_bo(x[cond_inside])})
        data_inside["es"] = "bo"
        data = pd.concat([data, data_inside], axis=0)
    if nb_above > 0:
        data_above = pd.DataFrame({'z': x[cond_above],
                                   'tau': x_tau[cond_above],
                                   'beta': x_beta[cond_above],
                                   'u': x_u[cond_above],
                                   'b1': 1.0 * np.ones(nb_above),
                                   'bo': 1.0 * np.ones(nb_above)})
        data_above["es"] = "pspl"
        data = pd.concat([data, data_above], axis=0)
    if nb_below > 0:
        data_below = pd.DataFrame({'z': x[cond_below],
                                   'tau': x_tau[cond_below],
                                   'beta': x_beta[cond_below],
                                   'u': x_u[cond_below],
                                   'b1': 0.0 * np.ones(nb_below),
                                   'bo': 0.0 * np.ones(nb_below)})
        data_below["es"] = "bo"
        data = pd.concat([data, data_below], axis=0)

    data = data.sort_values(['tau'])
    return data
#
# ========== Amplification ==========
def amplification_limb(t0, tE, u0, rho, gamma, t, piEN, piEE, DsN, DsE,
        filename_bo, filename_b1, flag_plot=0, path_output=".", t_pattern='n'):

    # Load tabulated Bo, interpolation and extrapolation
    file_Bo = open(filename_bo, 'r')
    file_Bo_p = pickle.load(file_Bo)
    z_bo = file_Bo_p.T[0]
    bo = file_Bo_p.T[1]
    file_Bo.close()
    f_bo = interp1d(z_bo, bo, kind='linear')

    file_B1 = open(filename_b1, 'r')
    file_B1_p = pickle.load(file_B1)
    z_b1 = file_B1_p.T[0]
    b1 = file_B1_p.T[1]
    file_B1.close()
    f_b1 = interp1d(z_b1, b1, kind='linear')

    # Amplification
    tau = (t-t0)/tE + piEN * DsN + piEE * DsE
    beta = u0 + piEN * DsE - piEE * DsN

    tau = (t-t0)/tE
    # beta = u0 * np.ones(t.shape[0])
    u = np.sqrt(tau**2 + beta**2)

    z_curr = u/rho
    data = f_bo_extrapol(z_curr, np.min(z_bo), np.max(z_bo), tau, beta, u, f_bo, f_b1)
    data["amp"] = (data.bo - gamma*data.b1) * (data.u**2 + 2)/(data.u * np.sqrt(data.u**2 + 4))

    togive = np.array([])
    for taucurr in tau:
        togive = np.append(togive, data[data["tau"]==taucurr]['amp'])

    return togive





















# -*-coding:Utf-8 -*
# ====================================================================

# ====================================================================
#   Packages
# ====================================================================
import sys
import numpy as np
from muLAn.models.vbb.vbb import vbbmagU
# ====================================================================
#   Functions
# ====================================================================
def magnifcalc(t, param, Ds=None, tb=None):
    """Return the VBB method finite-source uniform magnification."""
### Get parameters
    t0 = param['t0']
    u0 = param['u0']
    tE = param['tE']
    rho = param['rho']
    q = param['q']
    piEN = param['piEN']
    piEE = param['piEE']
    alpha0 = param['alpha']
    s0 = param['s']
    dalpha = param['dadt']
    ds = param['dsdt']
### Lens orbital motion
    alpha, s = lens_rotation(alpha0, s0, dalpha, ds, t, tb)
### Parallax
    DsN = Ds['N']
    DsE = Ds['E']
    tau = (t-t0)/tE + piEN * DsN + piEE * DsE
    beta = u0 + piEN * DsE - piEE * DsN
    x, y = binrot(alpha, tau, beta)
### Conversion secondary body left -> right
    x = -x
### Compute magnification
    accuracy = 1.e-2 # Absolute mag accuracy (mag+/-accuracy)
    return np.array([vbbmagU(s[i], q, rho, x[i], y[i], accuracy) for i in range(len(x))])
# --------------------------------------------------------------------
def binrot(theta, x_old, y_old):
    """Rotation by an angle alpha.

    :param theta: float, angle in radians.
    :param x_old: numpy array, x coodinate in the old frame.
    :param y_old: numpy array, y coodinate in the old frame.
    :return x_new: numpy array, x coodinate in the new frame.
    :return y_new: numpy array, y coodinate in the new frame.
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_new = x_old * cos_theta - y_old * sin_theta
    y_new = x_old * sin_theta + y_old * cos_theta
    return x_new, y_new
# --------------------------------------------------------------------
def lens_rotation(alpha0, s0, dalpha, ds, t, tb):
    """Compute the angle alpha and projected separation s for each
    time step due to the lens orbital motion.

    :param alpha0: angle alpha at date tb.
    :param s0: projected separation at date tb.
    :param dalpha: float, angular velocity at date tb
        (radians.year^-1).
    :param ds: change rate of separation (year^-1).
    :param t: list of dates.
    :param tb: time reference for linear development.
    :type alpha0: float
    :type s0: float
    :type dalpha: float
    :type ds: float
    :type t: numpy array
    :type tb: float
    :return: unpacked list of actual alpha and s values at each date.
    :rtype: numpy array, numpy array
    """
    Cte_yr_d = 365.25  # Julian year in days
    alpha = alpha0 - (t - tb) * dalpha / Cte_yr_d
    s = s0 + (t-tb) * ds / Cte_yr_d
    return alpha, s


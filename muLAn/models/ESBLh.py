# -*-coding:Utf-8 -*
# ====================================================================

# ====================================================================
# Standard packages
# ====================================================================
import sys
import numpy as np
# ====================================================================
# Non-standard packages
# ====================================================================
import esbltaylor.esbltaylor as esbltaylorbuilt
# ====================================================================
#   Functions
# ====================================================================
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
# --------------------------------------------------------------------
def magnifcalc(t, param, Ds=None, tb=None, **kwargs_method):
    # Compute the amplification
    kwargs = dict()
    kwargs.update(kwargs_method)
    kwargs.update({'params': param})
    kwargs.update({'dates': t})
    kwargs.update({'tb': tb})
    kwargs.update({'Ds': Ds})
    kwargs.update({'degree': 4})
    amp, flag = magnifcalc_wrap(**kwargs)
    return amp
# --------------------------------------------------------------------
def magnifcalc_wrap(**kwargs):

    try:
        params = kwargs['params']
    except KeyError:
        chat = "No parameters received in magnifcalc(...) from test_import."
        sys.exit(chat)
    try:
        t = kwargs['dates']
    except KeyError:
        chat = "No dates received in magnifcalc(...) from test_import."
        sys.exit(chat)
    try:
        tb = kwargs['tb']
    except KeyError:
        tb = params['t0']
    try:
        Ds = kwargs['Ds']
    except KeyError:
        Ds = dict({'N' : np.zeros(len(t)), 'E' : np.zeros(len(t))})
    try:
        degree = kwargs['degree']
    except KeyError:
        degree = 0
    try:
        err = float(kwargs['TriggerNextMethod'.lower()])
    except KeyError:
        err = 1e-3
    try:
        ray_sigma = float(kwargs['PrecisionGoalRayshooting'.lower()])
    except KeyError:
        ray_sigma = 1e-2
    try:
        ray_rect_pix = int(kwargs['LocalMarginRayshooting'.lower()])
    except KeyError:
        ray_rect_pix = 1

    t0 = params['t0']
    u0 = params['u0']
    tE = params['tE']
    rho = params['rho']
    gamma = params['gamma']
    q = params['q']
    piEN = params['piEN']
    piEE = params['piEE']
    alpha0 = params['alpha']
    s0 = params['s']
    dalpha = params['dadt']
    ds = params['dsdt']

    # Correction of the separation/angle due to lens orbital motion
    alpha, s = lens_rotation(alpha0, s0, dalpha, ds, t, tb)

    # Correction of the trajectory due to parallax
    DsN = Ds['N']
    DsE = Ds['E']
    tau = (t-t0)/tE + piEN * DsN + piEE * DsE
    beta = u0 + piEN * DsE - piEE * DsN
    x, y = binrot(alpha, tau, beta)

    # Center of mass to Cassan (2008)
    GL1 = s * q / (1 + q)
    x = x - GL1

    # Compute magnification using PSBL
    list = [[q, rho, gamma, degree, err, ray_sigma, ray_rect_pix], s.tolist(), x.tolist(), y.tolist()]
    magnif = esbltaylorbuilt.magnifcalc(list)

    return np.array(magnif[0]), np.array(magnif[1])


# ====================================================================
# Main function
# ====================================================================
if __name__=="__main__":
    print 'Please use with muLAn.'

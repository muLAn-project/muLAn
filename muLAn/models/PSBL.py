# -*-coding:Utf-8 -*
# ====================================================================

# ====================================================================
#   Packages
# ====================================================================
import sys
import numpy as np
# ====================================================================
#   Functions
# ====================================================================
def magnifcalc(t, param, Ds=None, tb=None):
    """Return the quadrupolar approximation of the magnification."""
### Get parameters
    t0 = param['t0']
    u0 = param['u0']
    tE = param['tE']
    q = param['q']
    piEN = param['piEN']
    piEE = param['piEE']
    alpha0 = param['alpha']
    s0 = param['s']
    dalpha = param['dalpha']
    ds = param['ds']
### Lens orbital motion
    alpha, s = lens_rotation(alpha0, s0, dalpha, ds, t, tb)
### Parallax
    DsN = Ds['N']
    DsE = Ds['E']
    tau = (t-t0)/tE + piEN * DsN + piEE * DsE
    beta = u0 + piEN * DsE - piEE * DsN
    x, y = binrot(alpha, tau, beta)
### Conversion center of mass to Cassan (2008)
    x = x - s*q/(1.+q)
### Compute magnification
    zeta0 = x + 1j*y
    return [monopole(s[i], q, zeta0[i]) for i in range(len(x))]
# --------------------------------------------------------------------
def monopole(s, q, zeta0):
    z0 = solve_lens_poly(s, q, zeta0) # convention Cassan (2008)
    W1 = 1./(1.+q)*(1./z0+q/(z0+s))
    z0 = z0[np.abs(z0-W1.conjugate()-zeta0)<0.000001]
    W2 = -1./(1.+q)*(1./z0**2+q/(z0+s)**2)
    mu0 = 1./(1.-np.abs(W2)**2)
    A0 = np.sum(np.abs(mu0))
    return A0
# --------------------------------------------------------------------
def solve_lens_poly(s,q,zeta):
    """Solve binary lens equation [convention Cassan (2008)]."""
    coefs = [(1+q)**2*(s+zeta.conjugate())*zeta.conjugate(),(1+q)*(s*(q-abs(zeta)**2*(1+q))+(1+q)*((1+2*s**2)-abs(zeta)**2+2*s*zeta.conjugate())*zeta.conjugate()),(1+q)*(s**2*q-s*(1+q)*zeta+(2*s+s**3*(1+q)+s**2*(1+q)*zeta.conjugate())*zeta.conjugate()-2*abs(zeta)**2*(1+q)*(1+s**2+s*zeta.conjugate())),-(1+q)*(s*q+s**2*(q-1)*zeta.conjugate()+(1+q+s**2*(2+q))*zeta+abs(zeta)**2*(2*s*(2+q)+s**2*(1+q)*(s+zeta.conjugate()))),-s*(1+q)*((2+s**2)*zeta+2*s*abs(zeta)**2)-s**2*q,-s**2*zeta]
    return np.roots(coefs)
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


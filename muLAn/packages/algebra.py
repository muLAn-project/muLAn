# -*-coding:Utf-8 -*

# Packages
# ========
import ConfigParser as cp
import copy
import glob
import muLAn
import muLAn.packages.general_tools as gtools
import numpy as np
import os
import pandas as pd
from sklearn import linear_model
import sys

# Functions
# =========
def test():
    print "Je te vois"

def fsfbwsig(time_serie, cond, blending=True):
    """
    Compute the source and blend flux using a linear fit with errors
    in flux.

    Parameters
    ----------
    time_serie : dict
        A dictionnary with the data flux (`flux`) and errors
        ('err_flux'), the corresponding magnification from the model
        (`amp`).

    cond : array-like
        Array of booleans used as a mask for ``time_serie``.

    blending : bool, default True
        Wheter or not fitting with a blending.

    Returns
    -------
    fs : float
        A float which is the source flux.
    fb : float
        A float which is the blend flux.
    """

    x = np.atleast_2d(time_serie['amp'][cond]).T
    y = np.atleast_2d(time_serie['flux'][cond]).T
    sig = np.atleast_2d(time_serie['err_flux'][cond]).T
    x2 = np.power(x, 2)
    sig2 = np.power(sig, 2)

    if blending:
        try:
            s = np.sum(1.0 / sig2)
            sx = np.sum(x / sig2)
            sy = np.sum(y / sig2)
            sxx = np.sum(x2 / sig2)
            sxy = np.sum(x * y / sig2)
            den = s * sxx - sx**2
            fs = (s * sxy - sx * sy) / den
            fb = (sxx * sy - sx * sxy) / den
        except ZeroDivisionError:
            fs = np.inf
            fb = np.inf
    else:
        fs, fb = fsfb(time_serie, cond, blending=False)

    if (np.abs(fs) == np.inf) | (np.abs(fb) == np.inf):
        fs, fb = fsfb(time_serie, cond, blending=False)
        print fs, fb

    return fs, fb


def fsfb(time_serie, cond, blending=True):

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


# Direct execution
# ================
if (__name__ == "__main__"):
    pass


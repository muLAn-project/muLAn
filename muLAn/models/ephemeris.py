# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
#
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

# filename = full_path + '../' + '.pythonexternallibpath'
# file = open(filename, 'r')
# for line in file:
#     path_lib_ext=line
# file.close()
# if path_lib_ext != 'None':
#     sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#   Packages
# ----------------------------------------------------------------------
import numpy as np
import configparser as cp
import glob
from PyAstronomy import pyasl
from astropy.time import Time
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord

# ----------------------------------------------------------------------
#   Functions
# ----------------------------------------------------------------------
def dsndse(l, b, cfgsetup, obs_properties, observatories):

    # Observer position
    path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

    if len(obs_properties['loc']) > 1:
        name1 = obs_properties['loc'][np.where(np.array(
            [obs == cfgsetup.get('Observatories', 'Reference').lower()
             for obs in observatories]) == True)[0][0]]
        name1 = glob.glob(path + name1 + '.*')[0]
    else:
        name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]


    name2 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]
    try:
        sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
            Ds(name1, name2, l, b,
                    cfgsetup.getfloat('Modelling', 'tp'), cfgsetup)
        if name1 != name2:
            DsN = DsNs
            DsE = DsEs
        else:
            DsN = DsNe
            DsE = DsEe
    except:
        DsN = None
        DsE = None

    return DsN, DsE


def Ds(ephearth, ephsat, l, b, tp, cfgsetup):
    ''' Compute Delta(s) of the Sun projected into the plane of the sky (T,E,N).
        Usage:
        - Annual parallax: Ds("Hori-Earth_2012-17.dat","Hori-Earth_2012-17.dat")
        - Spitzer space parallax: Ds("Hori-Earth_2012-17.dat","Hori-Spitzer_2012-17.dat")
        - K2 space parallax: Ds("Hori-Earth_2012-17.dat","Hori-K2_2012-17.dat")
        Returns:
        - (s,Ds)(T,E,N) as Earth-to-Sun or Sat-to-Earth vectors '''
    hjd,s3xe,s3ye,s3ze,v3xe,v3ye,v3ze = np.loadtxt(ephearth,usecols=(0,5,6,7,8,9,10),dtype='Float64',unpack=True)

    # t = hjd-2450000.
    # Time conversion: TDB->TCG->HJD
    c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), dec=cfgsetup.get('EventDescription', 'DEC'), frame='icrs')
    t = hjd-2400000.0
    ttest = t-50000.0
    flag_clem = 0
    if flag_clem:
        t = np.array([pyasl.helio_jd(tc, c_icrs.ra.degree, c_icrs.dec.degree) for tc in Time(t, format='mjd', scale='tdb').tcg.value])-50000.0
    else:
        t = np.array([pyasl.helio_jd(tc, l,b) for tc in Time(t, format='mjd', scale='tdb').tcg.value])-50000.0
    # print np.mean(t-ttest), np.min(t-ttest), np.max(t-ttest)

    s3xp = interp1d(t,s3xe,kind='linear')(tp)
    s3yp = interp1d(t,s3ye,kind='linear')(tp)
    s3zp = interp1d(t,s3ze,kind='linear')(tp)
    v3xp = interp1d(t,v3xe,kind='linear')(tp)
    v3yp = interp1d(t,v3ye,kind='linear')(tp)
    v3zp = interp1d(t,v3ze,kind='linear')(tp)
    Ds3x = s3xe - s3xp - np.multiply(t-tp,v3xp)
    Ds3y = s3ye - s3yp - np.multiply(t-tp,v3yp)
    Ds3z = s3ze - s3zp - np.multiply(t-tp,v3zp)
    refDsT,refDsE,refDsN = elcip2ciel(l,b,Ds3x,Ds3y,Ds3z)
    DsTe = interp1d(t,refDsT,kind='linear')
    DsEe = interp1d(t,refDsE,kind='linear')
    DsNe = interp1d(t,refDsN,kind='linear')
    refsT,refsE,refsN = elcip2ciel(l,b,s3xe,s3ye,s3ze)
    sTe = interp1d(t,refsT,kind='linear')
    sEe = interp1d(t,refsE,kind='linear')
    sNe = interp1d(t,refsN,kind='linear')

    hjd, s3x, s3y, s3z = np.loadtxt(ephsat, usecols=(0,5,6,7), dtype='Float64', unpack=True)

    # t = hjd-2450000.
    # Time conversion: TDB->TCG->HJD
    t = hjd-2400000.0
    if flag_clem:
        t = np.array([pyasl.helio_jd(tc, c_icrs.ra.degree, c_icrs.dec.degree) for tc in Time(t, format='mjd', scale='tdb').tcg.value])-50000.0
    else:
        t = np.array([pyasl.helio_jd(tc, l,b) for tc in Time(t, format='mjd', scale='tdb').tcg.value])-50000.0

    s3xs = s3xe - s3x
    s3ys = s3ye - s3y
    s3zs = s3ze - s3z
    Ds3x = s3xs - s3xp - np.multiply(t-tp,v3xp)
    Ds3y = s3ys - s3yp - np.multiply(t-tp,v3yp)
    Ds3z = s3zs - s3zp - np.multiply(t-tp,v3zp)
    refDsT,refDsE,refDsN = elcip2ciel(l,b,Ds3x,Ds3y,Ds3z)
    DsTs = interp1d(t,refDsT,kind='linear')
    DsEs = interp1d(t,refDsE,kind='linear')
    DsNs = interp1d(t,refDsN,kind='linear')
    refsT,refsE,refsN = elcip2ciel(l,b,s3xs,s3ys,s3zs)
    sTs = interp1d(t,refsT,kind='linear')
    sEs = interp1d(t,refsE,kind='linear')
    sNs = interp1d(t,refsN,kind='linear')

    return sTe,sEe,sNe,DsTe,DsEe,DsNe,sTs,sEs,sNs,DsTs,DsEs,DsNs

def elcip2ciel(l,b,Xe,Ye,Ze):
    ''' Convert ecliptic coordinates (Xe,Ye,Ze) to plane of the sky (T,E,N) coordinates. '''
    l = l*np.pi/180.
    b = b*np.pi/180.
    cosl = np.cos(l)
    cosb = np.cos(b)
    sinl = np.sin(l)
    sinb = np.sin(b)
    T = Xe*cosb*cosl + Ye*cosb*sinl + Ze*sinb # T = Target
    E = -Xe*sinl + Ye*cosl # E = East
    N = -Xe*sinb*cosl -Ye*sinb*sinl + Ze*cosb # N = North
    return T,E,N

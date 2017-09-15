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

def magnifcalc(timeserie, physical_properties, Ds=0, tb=None):

    u0  = physical_properties['u0']
    tE = physical_properties['tE']
    t0 = physical_properties['t0']
    piEN = physical_properties['piEN']
    piEE = physical_properties['piEE']

    DsN = Ds['N']
    DsE = Ds['E']

    tau = (timeserie-t0)/tE + piEN * DsN + piEE * DsE
    beta = u0 + piEN * DsE - piEE * DsN

    return magnification(tau, beta)

# ========= Magnification calculation =========
def magnification(tau, beta):
    u = np.sqrt(tau**2 + beta**2)
    A = (u**2+2)/(u*np.sqrt(u**2+4))
    return A

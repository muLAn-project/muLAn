# -*-coding:Utf-8 -*
# ====================================================================
# Packages
# ====================================================================
import configparser as cp
import copy
import glob
import muLAn
import muLAn.packages.general_tools as gtools
import muLAn.packages.algebra as algebra
import numpy as np
import os
import pandas as pd
import sys
import tables

class Instrument:
    """Class to store properties of each instrument

    Args:
        properties (list): list with shape (1, 2), where properties[0] is
            the ID of the observatory (str), and properties[1] is a string
            describing the properties of the observatory, as follow,
            "Label, HTML color, Passband[=Gamma], Type, Location[, list of int".
            `Label` is the name of the instrument; the color should be written
            without the `#`, Gamma is the linear limb darkening coefficient
            (default: 0.0), `Type` is `Magnitude` or `Flux` depending of the
            input data files used, `Location` is a name referring to the
            position of the instrument (a file `Location.dat` with JPL ephemeris
            is required), and the optional list of int correspond to the data ID
            to remove before the fit and the plots.

            Examples for properties[1]:

            OGLE-I, 000000, I=0.64, Magnitude, Earth
            OGLE-V, 000000, V=0.5, Magnitude, Earth, 11, 60, 68, 73, 78, 121, 125, 128, 135
            OGLE Passband I, 000000, I, Magnitude, Earth, 11, 60-68

            In the third example, the data points from 60 to 68 (included) will
            be removed. A file Earth.dat should be in the Data/ directory.
        
    Attributes:
        id (str): instrument ID.
        label (str): instrument label.
        color (str): HTML color with #.
        passband (str): passband name.
        gamma (float): linear limb-darkening coefficient Gamma.
        type (str): wether the input data are in flux or magnitude units.
        location (str): key word corresponding to the ephemeris file.
        reject (:obj:`numpy.array`): array of the data ID to remove.

    """

    def __init__(self, properties):

        self._extract_properties(properties)


    def _extract_properties(self, prop):

        properties = dict()
        properties['id'] = prop[0]
        props = prop[1].split(',')
        n_opts = len(props)
        if n_opts > 4:
            keywd = 'label color band type location'.split(' ')
            for i in range(5):
                properties.update({keywd[i]: props[i].strip()})

            if n_opts > 5:
                properties.update({'reject': props[5:]})
                self.reject = self._rejection_list(properties['reject'])
        else:
            txt = 'Syntax error or not enough properties provided for an instrument.'
            sys.exit(txt)


        self.id = properties['id']
        self.label = properties['label']
        self.color = '#{:s}'.format(properties['color'])
        band = self._extract_band(properties['band'])
        self.passband = band[0]
        self.gamma = band[1]
        self.type = properties['type']
        self.location = properties['location']


    def _rejection_list(self, string):

        string = [a.strip() for a in string]
        nb = len(string)
        to_reject = np.array([], dtype=np.int)
        for i in range(nb):
            substring = string[i].split('-')
            if len(substring) == 1:
                to_reject = np.append(to_reject, int(substring[0]))
            elif len(substring) == 2:
                a = int(substring[0].strip())
                b = int(substring[1].strip())
                n = b - a + 1
                ids = np.linspace(a, b, n, dtype=np.int)
                to_reject = np.append(to_reject, ids) 

        return to_reject


    def _extract_band(self, string):

        string = string.split('=')
        string = [a.strip() for a in string]
        if len(string) == 2:
            return string[0].strip(), float(string[1])
        else:
            return string[0].strip(), 0.0

class InstrumentsList(Instrument):
    """Class to store a list of instruments (observatories).

    Args:
        input (str): file that defines instrument properties.
    
    Attributes:
        to_dict(dict): dictionary with keys corresponding to each instrument
            ID, and values corresponding to a :obj:`muLAn.instruments.Instrument`
            object.

    """

    def __init__(self, input):

        if isinstance(input, str):
            self.file = input
            self._load_from_file(input)

        self._create_instruments_list()

    def _load_from_file(self, fname):

        cfgobs = cp.ConfigParser()
        cfgobs.read(fname)
        self.parser = cfgobs

    def _create_instruments_list(self):

        item = self.parser.items('ObservatoriesDetails')
        n = len(item)
        instruments_list = dict()
        for i in range(n):
            tmp = Instrument(list(item[i])) 
            instruments_list.update({tmp.id: tmp})
            setattr(self, tmp.id, tmp)

        self._instruments_list = instruments_list
        self.len = n

    def to_dict(self):

        return self._instruments_list

    def prop(self, val):

        return self._instruments_list[val[0]]



if (__name__ == "__main__"):
    pass


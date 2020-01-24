# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# Call the modules asked by the user
# ----------------------------------------------------------------------
# Packages
# ----------------------------------------------------------------------
from astropy.coordinates import SkyCoord
import configparser as cp
import glob
import importlib
from muLAn.data import Data
import muLAn.iotools as iotools
from muLAn.iotools import LensModel
import muLAn.models as mulanmodels
import muLAn.models.ephemeris as ephemeris
from muLAn.packages.general_tools import *
from muLAn.instruments import InstrumentsList
import muLAn.packages.sortmodels as mulansort
import muLAn.plottypes as mulanplots
import numpy as np
import pandas as pd
import os
import pickle
import scipy
import sys

# ====================================================================
# Fonctions
# ====================================================================
def run_sequence(path_event, options):

    # Load configuration files
    cfgsetup = cp.ConfigParser()
    cfgsetup.read([path_event + 'setup.ini', path_event + 'advancedsetup.ini'])

    text = "Load parameter files..."
    communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True)

    # Deprecated:
    cfgobs = cp.ConfigParser()
    cfgobs.read(path_event + 'observatories.ini')
    # Replaced by:
    fname = "{:s}/observatories.ini".format(path_event)
    instruments = InstrumentsList(fname)

    # Add the path to the configuration
    cfgsetup.set('FullPaths', 'Event', path_event)

    # Check the paths
    if cfgsetup.get('FullPaths', 'Code').replace(" ", "") != "":
        if cfgsetup.get('FullPaths', 'Code')[-1] != '/':
            cfgsetup.set('FullPaths', 'Code', cfgsetup.get('FullPaths', 'Code') + '/')
        if cfgsetup.get('FullPaths', 'Event')[-1] != '/':
            cfgsetup.set('FullPaths', 'Event', cfgsetup.get('FullPaths', 'Event') + '/')
    if cfgsetup.get('RelativePaths', 'Data')[-1] != '/':
        cfgsetup.set('RelativePaths', 'Data', cfgsetup.get('RelativePaths', 'Data') + '/')
    if cfgsetup.get('RelativePaths', 'Plots')[-1] != '/':
        cfgsetup.set('RelativePaths', 'Plots', cfgsetup.get('RelativePaths', 'Plots') + '/')
    if cfgsetup.get('RelativePaths', 'Chains')[-1] != '/':
        cfgsetup.set('RelativePaths', 'Chains', cfgsetup.get('RelativePaths', 'Chains') + '/')
    if cfgsetup.get('RelativePaths', 'Outputs')[-1] != '/':
        cfgsetup.set('RelativePaths', 'Outputs', cfgsetup.get('RelativePaths', 'Outputs') + '/')
    if cfgsetup.get('RelativePaths', 'Archives')[-1] != '/':
        cfgsetup.set('RelativePaths', 'Archives', cfgsetup.get('RelativePaths', 'Archives') + '/')
    if cfgsetup.get('RelativePaths', 'ModelsHistory')[-1] != '/':
        cfgsetup.set('RelativePaths', 'ModelsHistory', cfgsetup.get('RelativePaths', 'ModelsHistory') + '/')

    if not cfgsetup.has_section('Optimization'):
        cfgsetup.add_section('Optimization')
    if not cfgsetup.has_option('Optimization', 'UseBinaryFiles'):
        cfgsetup.set('Optimization', 'UseBinaryFiles', 'False')

    if not cfgsetup.has_section('Modelling'):
        cfgsetup.add_section('Modelling')
    if not cfgsetup.has_option('Modelling', 'IncludeBlending'):
        cfgsetup.set('Modelling', 'IncludeBlending', 'True')

    # Take into account manual options
    if options['plot'] != None:
        cfgsetup.set('Controls', 'Modes', 'Plot')
        cfgsetup.set('Plotting', 'Models', options['plot'])
    if options['fit']:
        cfgsetup.set('Controls', 'Modes', 'Fit')
    cond = (options['fit']) and (options['plot'] != None)
    if cond:
        cfgsetup.set('Controls', 'Modes', 'Fit, Plot')
        cfgsetup.set('Plotting', 'Models', options['plot'])
    if options['archive'] != None:
        cfgsetup.set('Controls', 'Archive', options['archive'])
    if options['ncores'] != None:
        if options['ncores'] > 0:
            cfgsetup.set('FitSetupDMCMC', 'Threads', '{:d}'.format(options['ncores']))
    if options['nchains'] != None:
        if options['nchains'] > 0:
            cfgsetup.set('FitSetupDMCMC', 'Chains', '{:d}'.format(options['nchains']))
    if options['resume'] != None:
        if options['resume']: cfgsetup.set('FitSetupDMCMC', 'Resume', 'True')
    if options['verbose'] != None:
        if options['verbose'] > -1:
            cfgsetup.set('Modelling', 'Verbose', '{:d}'.format(options['verbose']))
    if options['optimize'] != None:
        if options['optimize']: cfgsetup.set('Controls', 'Optimize', 'True')
        else: cfgsetup.set('Controls', 'Optimize', 'False')

    # Controls
    modes = [a.lower() for a in unpack_options(cfgsetup, 'Controls', 'Modes')]
    cfgsetup.set('Modelling', 'Fit', str(any(x == 'fit' for x in modes)))
    cfgsetup.set('Plotting', 'flag', str(any(x == 'plot' for x in modes)))

    # Verbose
    text = "Parameter files have been read and checked."
    communicate(cfgsetup, 4, text, opts=[printoption.good], prefix=False, newline=False)

    # Check directories and create missing ones
    paths_to_check = ['Data', 'Plots', 'Outputs', 'Chains', 'Archives']
    for i in range(len(paths_to_check)):
        if not os.path.exists(cfgsetup.get('RelativePaths', paths_to_check[i])):
            os.makedirs(cfgsetup.get('RelativePaths', paths_to_check[i]))

    # Test the archive
    if cfgsetup.getboolean('Modelling', 'Fit') & (cfgsetup.getboolean('FitSetupDMCMC', 'Resume')==False)\
            & (cfgsetup.getint('Modelling', 'Verbose')>0):
        filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Archives')
        name = prompt(filename, cfgsetup.get('Controls', 'Archive'))
        cfgsetup.set('Controls', 'Archive', name)

    # Print the details
    text = "Event details"
    communicate(cfgsetup, 3, text, opts=[printoption.level1], prefix=False, newline=True)

    text = cfgsetup.get("EventDescription", "Name")
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

    text = "Original equatorial coordinates:    {:s} {:s}".format(cfgsetup.get("EventDescription", "RA"), cfgsetup.get("EventDescription", "DEC"))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

    c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), dec=cfgsetup.get('EventDescription', 'DEC'), frame='icrs')
    text = "Equatorial coordinates [deg]:       {:.6f}   {:.6f}".format(c_icrs.ra.degree, c_icrs.dec.degree)
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
    b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree
    text = "Ecliptic coordinates [deg]:         {:10.6f}   {:10.6f}".format(l, b)
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    l = c_icrs.transform_to('galactic').l.degree
    b = c_icrs.transform_to('galactic').b.degree
    text = "Galactic coordinates [deg]:         {:10.6f}   {:10.6f}".format(l, b)
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    text = "Modelling options"
    communicate(cfgsetup, 3, text, opts=[printoption.level1], prefix=False, newline=True)

    models_temp = np.array(cfgsetup.options("Modelling"))
    cond = np.where(np.array([a.find("models_") for a in cfgsetup.options("Modelling")])==0)
    models_temp = np.atleast_1d(models_temp[cond])
    try:
        a = models_temp[0]
    except:
        raise KeyError("No models to load.")
    models_temp2 = np.atleast_1d([unpack_options(cfgsetup, "Modelling", a) for a in models_temp])
    for i in range(len(models_temp)):
        text = "Models for {:s} observations:".format(models_temp[i].split("_")[1].upper())
        communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)
        for j in range(len(models_temp2[i])):
            a = np.array(models_temp2[i][j].split('/'))
            text = "  {:s} {:.6f} --> {:.6f}".format(a[1], float(a[0].split("-")[0]), float(a[0].split("-")[1]))
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    text = "Minimisation: {:s}".format(cfgsetup.get("Modelling", "Method").upper())
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

    text = "List of parameters"
    communicate(cfgsetup, 3, text, opts=[printoption.level1], prefix=False, newline=True)

    counter = 0
    text = "q      {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "q")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "s      {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "s")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "tE     {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "tE")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "rho    {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "rho")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "piEN   {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "piEN")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "piEE   {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "piEE")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1

    text = "t0     {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "t0")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "u0     {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "u0")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "alpha  {:s}".format(print_params(unpack_options(cfgsetup, "Modelling", "alpha")))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "tp     {:s}".format(cfgsetup.get("Modelling", "tp"))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    if text.find("within") > -1: counter = counter + 1
    text = "tb   {:s}".format(cfgsetup.get("Modelling", "tb"))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    text = "Fit parameters"
    communicate(cfgsetup, 3, text, opts=[printoption.level1], prefix=False, newline=True)
    text = "Markov Chain Monte Carlo controls:"
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

    if (cfgsetup.getint("FitSetupDMCMC", "Chains")%2) | (cfgsetup.getint("FitSetupDMCMC", "Chains") < 2.0*counter):
        text = "\n\033[1m\033[91mThe chains number must be even and at least twice the number of fitted\nparameters. muLAn killed.\033[0m"
        sys.exit(text)
    else:
        text = "  {:s} chains in parallel".format(cfgsetup.get("FitSetupDMCMC", "Chains"))
        communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    text = "  Length of each chain: {:s}".format(cfgsetup.get("FitSetupDMCMC", "ChainLength"))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
    text = "  Maximum number of threads: {:s}".format(cfgsetup.get("FitSetupDMCMC", "Threads"))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    text = "Plot options:"
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

    plots_temp = np.atleast_1d(cfgsetup.get("Plotting", "Type"))
    options_temp = np.atleast_1d(cfgsetup.get("Plotting", "Options"))
    text = "  "
    [communicate(cfgsetup, 3, text + plots_temp[i] + ": " + options_temp[i], opts=False, prefix=False, newline=False, tab=True) for i in range(len(plots_temp))]


    # ----------------------------------------------------------------------
    #   Data
    # ----------------------------------------------------------------------
    if cfgsetup.getboolean('Plotting', 'Data') \
            | cfgsetup.getboolean('Modelling', 'Fit'):

        # communicate(cfgsetup, 1, "\n")
        # text = " Data loading & fit preparation "
        # communicate(cfgsetup, 1, text, opts=[printoption.reverse])

        text = "Load data"
        communicate(cfgsetup, 3, text, opts=[printoption.level0], prefix=True, newline=True)

        # Data file names
        data2find = np.array(cfgobs.options('ObservatoriesDetails'))

        # Cross-match with data to use
        data2use = np.array(cfgsetup.options('Observatories'))
        data2use = np.array([data2find[i] for i in range(len(data2find)) if np.any(data2use == data2find[i])])

        # Reference for plotting
        filename = cfgsetup.get('FullPaths', 'Event') + 'setup.ini'
        file = open(filename, 'r')
        for line in file:
            a = line.split(':')[0]
            if a != '':
                a = a.strip().lower()
                if len(np.where(data2use == a)[0]) > 0:
                    cfgsetup.set('Observatories', 'Reference', a)
                    break

        # Cross-match with existing files
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

        data_filenames = glob.glob(path + '*')
        observatories = [a.split('/')[-1] for a in data_filenames]
        data_filenames = [data_filenames[i] for i in range(len(data_filenames)) if
                          np.any(data2use == observatories[i].rpartition('.')[0].lower())]
        observatories = [ob.rpartition('.')[0].lower() for ob in observatories if
                         np.any(data2use == ob.rpartition('.')[0].lower())]

        obs_properties = dict({'name': [], 'colour': [], 'filter': [], 'gamma': [],
                'loc': [], 'fluxoumag': [], 'exclude': [], 'key': []})
        # text = 'Data from:\n'
        for i in range(len(observatories)):
            table = [a.strip() for a in cfgobs.get('ObservatoriesDetails', observatories[i]).split(',')]

            obs_properties['name'].append(table[0])
            obs_properties['colour'].append(table[1])
            obs_properties['fluxoumag'].append(table[3])
            obs_properties['loc'].append(table[4])
            obs_properties['key'].append(observatories[i])

            try:
                obs_properties['filter'].append(table[2].split("=")[0])
                obs_properties['gamma'].append(float(table[2].split("=")[1]))
            except:
                obs_properties['filter'].append(table[2])
                obs_properties['gamma'].append([0.0])

            if (len(table[5:]) > 0):
                if (table[5] != ''):
                    obs_properties['exclude'].append(table[5:])
                else:
                    obs_properties['exclude'].append(None)
            else:
                obs_properties['exclude'].append(None)

        #     if i!=len(observatories)-1:
        #         text = text + '    ' + table[0] + '\n'
        #     if i==len(observatories)-1:
        #         text = text + '    ' + table[0]
        #
        # communicate(cfgsetup, 1, text)

        # Parameters of modelling
        model_params = dict()
        interpol_method = dict()
        for i in range(len(obs_properties['loc'])):
            name = 'Models_' + obs_properties['loc'][i]
            table = np.array([a.split('/')[1].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            model_params.update({name: table})

            table = np.array([a.split('/')[0].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            name2 = 'DateRanges_' + obs_properties['loc'][i]
            model_params.update({name2: table})

            for j in range(len(unpack_options(cfgsetup, 'Modelling', name))):
                a = unpack_options(cfgsetup, 'Modelling', name)[j]
                a = a.split('/')
                if len(a) == 3:
                    name3 = obs_properties['loc'][i] + '#' + a[1] + '#' + a[0]
                    t1 = float(a[0].split('-')[0].strip())
                    t2 = float(a[0].split('-')[1].strip())
                    # print(name3, t1, t2)
                    interpol_method.update({name3: [np.linspace(t1, t2, int(a[2])), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8')]})
                    # print(interpol_method)

        # Ephemeris
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

        if len(obs_properties['loc']) > 1:
            name1 = obs_properties['loc'][np.where(np.array(
                [obs == cfgsetup.get('Observatories', 'Reference').lower()
                 for obs in observatories]) == True)[0][0]]
            name1 = glob.glob(path + name1 + '.*')[0]
        else:
            name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]

        # Load data
        format = {'names': ('id', 'dates', 'magnitude', 'err_magn', 'seeing', 'background'), \
                  'formats': ('i8', 'f8', 'f8', 'f8', 'f8', 'f8')}

        # Event coordinates conversion
        c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), \
                          dec=cfgsetup.get('EventDescription', 'DEC'),
                          frame='icrs')
        l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
        b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

        # Use Pandas DataFrame to store data
        # ==================================
        filename_bin = "{:s}/all_data.bin".format(cfgsetup.get('RelativePaths', 'Data'))
        if (not cfgsetup.getboolean('Optimization', 'UseBinaryFiles'))\
                | (not os.path.exists(filename_bin)):

            data = pd.DataFrame()
            namecol = np.array('id dates magnitude err_magn_orig seeing background'.split())
            namecolf = np.array('id dates flux err_flux_orig seeing background'.split())
            coltype = 'i8 f8    f8        f8       f8    f8'.split()
            coltypet = np.array([np.dtype(a) for a in coltype])
            usecols = range(5)
            coltype = dict(zip(namecol, coltypet))
            coltypef = dict(zip(namecolf, coltypet))
            li = []
            model2load = np.array([])
            text_nb_data = ""

            # Identify the datasets providing flux
            list_flux = np.array([a for a in observatories 
                for i in range(len(obs_properties['key'])) 
                if ((a == obs_properties['key'][i]) & 
                    ((obs_properties['fluxoumag'][i]).lower() == "flux"))])

            for i in range(len(observatories)):
                # Load data
                flag_flux = any(list_flux == observatories[i])
                if flag_flux:
                    df = pd.read_csv(data_filenames[i], sep='\s+', names=namecolf,
                        usecols=usecols, dtype=coltypef, skiprows=1)
                else:
                    df = pd.read_csv(data_filenames[i], sep='\s+', names=namecol,
                        usecols=usecols, dtype=coltype, skiprows=1)

                # Add dataset name
                df['obs'] = observatories[i]
                nb_data_init = len(df)
                
                # Linear limb darkening coefficient Gamma
                df['gamma'] = obs_properties['gamma'][i]
                
                # Add location for ephemeris
                df['loc'] = obs_properties['loc'][i]
                
                # Rescale the error-bars [see Wyrzykowski et al. (2009)]
                gamma = float(unpack_options(cfgsetup, 'Observatories', observatories[i])[0][1:])
                epsilon = float(unpack_options(cfgsetup, 'Observatories', observatories[i])[1][:-1])
                if flag_flux:
                    df['err_flux'] = np.sqrt(np.power(gamma * df['err_flux_orig'],2)\
                            + np.power((np.log(10) * epsilon * df['flux']) / 2.5, 2))
                else:
                    df['err_magn'] = np.sqrt(np.power(gamma * df['err_magn_orig'], 2) + epsilon ** 2)

                # Select data only in specified time interval
                prop_temp = cfgsetup.get('Observatories', observatories[i]).split(',')

                if len(prop_temp) < 3:
                    txt = ("Format error in setup.ini [Observatories]"
                            "Name : (gamma, epsilon), time1-time2")
                    sys.exit(txt)
                prop_temp = prop_temp[2].replace(' ', '')
                if not prop_temp == '':
                    prop_temp = prop_temp.split('-')
                    try:
                        limits = [float(prop_temp[0]), float(prop_temp[1])]
                    except ValueError as err :
                        txt = "{:s}\n{:s}\n{:s}\n{:s} {:s}.".format(
                            "Format error in setup.ini",
                            "In [Observatories] you should write:",
                            "Name : (gamma, epsilon), time1-time2",
                            "No data will be removed from",
                            data_filenames[i])
                        print(err, txt)
                        continue

                    mask = (df['dates'] < limits[0]) | (df['dates'] >= limits[1])
                    df.drop(df[mask].index, inplace=True)

                # Remove outliers specified in Observatories.ini
                prop_temp = cfgobs.get('ObservatoriesDetails', observatories[i]).split(',')
                if len(prop_temp) > 5:
                    idx = np.array([], dtype='i8')
                    if prop_temp[5].strip() != '':
                        list_temp = [prop_temp[iii].strip() for iii in range(len(prop_temp)) if
                                     (iii > 4) & (prop_temp[iii].strip() != '')]
                        for a in list_temp:
                            toremove = np.array(a.split('-'), dtype='i8')
                            if len(toremove) == 1: 
                                idx = np.append(idx, toremove)
                            elif len(toremove) == 2: 
                                n = toremove[1] - toremove[0]
                                idx = np.append(idx, np.linspace(toremove[0], toremove[1], n+1, 
                                    endpoint=True, dtype='i8'))
                    
                        for j in range(idx.shape[0]):
                            mask = df['id'] == idx[j]
                            df.drop(df[mask].index, inplace=True)
                
                a = float(unpack_options(cfgsetup, "Observatories", observatories[i])[2].split("-")[0])
                b = float(unpack_options(cfgsetup, "Observatories", observatories[i])[2].split("-")[1])
                text = "Reading data within {:.6f} --> {:.6f} from {:s}".format(a, b, data_filenames[i].split("/")[-1])
                communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
                a = len(df)
                b = nb_data_init - a
                if b==0:
                    text_nb_data = text_nb_data + "    \033[1m{:6d}\033[0m data\n".format(a, b)
                else:
                    text_nb_data = text_nb_data + "    \033[1m{:6d}\033[0m data + \033[40m\033[97m\033[1m{:6d} data excluded\033[0m\n".format(a, b)
                if i==len(observatories)-1:
                    text_nb_data = text_nb_data + "  = \033[1m{:6d}\033[0m data in total".format(len(data))

                # Calculations from ephemeris for parallax
                name2 = glob.glob(path + obs_properties['loc'][i] + '.*')[0]
                try:
                    sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
                        ephemeris.Ds(name1, name2, l, b,
                                cfgsetup.getfloat('Modelling', 'tp'), cfgsetup)
                    if name1 != name2:
                        DsN = DsNs
                        DsE = DsEs
                    else:
                        DsN = DsNe
                        DsE = DsEe
                    df['DsN'] = DsN(df['dates'].values)
                    df['DsE'] = DsE(df['dates'].values)
                except:
                    DsN = None
                    DsE = None
                    df['DsN'] = 0
                    df['DsE'] = 0

                # Models
                name = 'Models_' + obs_properties['loc'][i]
                models_temp = model_params[name]
                name = 'DateRanges_' + obs_properties['loc'][i]
                dates_temp = model_params[name]

                df['model'] = 'PSPL'

                key = np.array([key for key in interpol_method])
                for j in range(len(models_temp)):
                    model2load = np.append(model2load, models_temp[j])
                    tmin = float((dates_temp[j]).split('-')[0].strip())
                    tmax = float((dates_temp[j]).split('-')[1].strip())

                    mask = (df['dates'] > tmin) & (df['dates'] <= tmax)
                    df.loc[mask, 'model'] = models_temp[j]

                if flag_flux:
                    try:
                        df['magnitude'] = 18.0 - 2.5 * np.log10(df['flux'])
                        df['err_magn_orig'] = np.abs(2.5 * df['err_flux_orig']\
                                / (df['flux'] * np.log(10)))
                        df['err_magn'] = np.sqrt(np.power(gamma * df['err_magn_orig'], 2) + epsilon ** 2)
                    except:
                        df['magnitude'] = 0.0
                        df['err_magn_orig'] = 0.0
                        df['err_magn'] = 0.0
                else:
                    # Compute flux from magnitudes
                    df['flux'] = np.power( 10, 0.4 * (18.0 - df['magnitude']) )
                    df['err_flux_orig'] = np.abs((np.log(10) / 2.5) *
                            df['err_magn_orig'] * df['flux'])
                    df['err_flux'] = np.abs((np.log(10) / 2.5) * df['err_magn'] * df['flux'])

                li.append(df)

            # Display data removed and included
            communicate(cfgsetup, 3, text_nb_data, opts=False, prefix=False, newline=True, tab=False)

            # Concatenate all the data
            data = pd.concat(li, axis=0, ignore_index=True)
            data.astype({'gamma': np.dtype('f8')})

            # Correct dates
            mask = data['dates'] > 2450000
            data.loc[mask, 'dates'] = data.loc[mask, 'dates'] - 2450000

            # Create columns for use in fit
            data['amp'] = -1
            data['fs'] = -999
            data['fb'] = -999
            li = ['amp', 'fs', 'fb']
            [data.astype({a: np.dtype('f8')}) for a in li]

            # For compatibility, convert DataFrame data into a dictionnary
            time_serie = data.to_dict('list')

            # Decide if method interpolation is used or not.
            key_list = np.array([key for key in interpol_method])
            if len(key_list) > 0:

                text = "Ask interpolation"
                communicate(cfgsetup, 3, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

                for i in range(len(key_list)):
                    loc = key_list[i].split('#')[0]

                    tmin = float(key_list[i].split('#')[2].split('-')[0])
                    tmax = float(key_list[i].split('#')[2].split('-')[1])

                    cond1 = (time_serie['dates'] <= tmax) & (time_serie['dates'] >= tmin) &\
                            (time_serie['loc'] == loc)

                    # print(loc, len(time_serie['model'][cond1]), len(interpol_method[key_list[i]][0]))

                    if len(time_serie['model'][cond1]) > len(interpol_method[key_list[i]][0]):

                        text = "Relevant for {:s} within {:s} --> {:s}\n".format(
                                loc,
                                key_list[i].split('#')[2].split('-')[0],
                                key_list[i].split('#')[2].split('-')[1])
                        text = text\
                                + "        {:d} observations / {:d} points asked: \033[1m\033[32mOK\033[0m".format(
                                len(time_serie['model'][cond1]),
                                len(interpol_method[key_list[i]][0]))
                        # communicate(cfgsetup, 1, text, opts=False)
                        communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

                        time_serie['interpol'][cond1] = key_list[i]

                        # Ephemeris
                        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

                        if len(obs_properties['loc']) > 1:
                            name1 = obs_properties['loc'][np.where(np.array(
                                [obs == cfgsetup.get('Observatories', 'Reference').lower()
                                 for obs in observatories]) == True)[0][0]]
                            name1 = glob.glob(path + name1 + '.*')[0]
                        else:
                            name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]

                        c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), \
                                          dec=cfgsetup.get('EventDescription', 'DEC'), frame='icrs')
                        # print(c_icrs.transform_to('barycentrictrueecliptic'))
                        l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
                        b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

                        name2 = glob.glob(path + loc + '.*')[0]
                        sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
                            ephemeris.Ds(name1, name2, l, b, cfgsetup.getfloat('Modelling', 'tp'), \
                                         cfgsetup)

                        if name1 != name2:
                            DsN = DsNs
                            DsE = DsEs
                        else:
                            DsN = DsNe
                            DsE = DsEe

                        interpol_method[key_list[i]][1] = DsN(interpol_method[key_list[i]][0])
                        interpol_method[key_list[i]][2] = DsE(interpol_method[key_list[i]][0])
                    else:
                        text = "Not relevant for {:s} within {:s} --> {:s}\n".format(
                                loc,
                                key_list[i].split('#')[2].split('-')[0],
                                key_list[i].split('#')[2].split('-')[1])
                        text = text\
                                + "        {:d} observations / {:d} points asked: \033[1m\033[31mIGNORED\033[0m".format(
                                len(time_serie['model'][cond1]),
                                len(interpol_method[key_list[i]][0]))
                        communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

                        del interpol_method[key_list[i]]

            filename = "{:s}/all_data.bin".format(cfgsetup.get('RelativePaths', 'Data'))
            outfile = open(filename,'wb')
            to_dump = [time_serie, model2load]
            pickle.dump(to_dump, outfile)
            outfile.close()

        # (Optimization) Load saved data from binary files
        elif cfgsetup.getboolean('Optimization', 'UseBinaryFiles')\
                & os.path.exists(filename_bin):

            text = "Data from the previous run, loaded from a binary file."
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
            text = "If you want to reload all the data, please use option UseBinaryFiles=False in advancedsetup.ini."
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

            filename = "{:s}/all_data.bin".format(cfgsetup.get('RelativePaths', 'Data'))
            infile = open(filename,'rb')
            list_input = pickle.load(infile)
            time_serie = list_input[0]
            model2load = list_input[1]
            infile.close()

        # Identify which quantities must be fit
        # -------------------------------------
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
            'dadt': np.array([a.strip() for a in cfgsetup.get('Modelling', 'dadt').split(',')]),\
            'dsdt': np.array([a.strip() for a in cfgsetup.get('Modelling', 'dsdt').split(',')])\
            }
        
        constp = dict()
        fitp = dict()
        result = np.array([])
        if params['t0'][0]!="fit":
            if params['t0'][0]=="gri":
                constp.update({'t0': node['t0']})
            else:
                constp.update({'t0': params['t0'][3].astype(np.float64)})
        else:
            fitp.update({'t0': params['t0'][3].astype(np.float64)})
            result = np.append(result, fitp['t0'])
        if params['u0'][0]!="fit":
            if params['u0'][0]=="gri":
                constp.update({'u0': node['u0']})
            else:
                constp.update({'u0': params['u0'][3].astype(np.float64)})
        else:
            fitp.update({'u0': params['u0'][3].astype(np.float64)})
            result = np.append(result, fitp['u0'])
        if params['tE'][0]!="fit":
            if params['tE'][0]=="gri":
                constp.update({'tE': node['tE']})
            else:
                constp.update({'tE': params['tE'][3].astype(np.float64)})
        else:
            fitp.update({'tE': params['tE'][3].astype(np.float64)})
            result = np.append(result, fitp['tE'])
        if params['rho'][0]!="fit":
            if params['rho'][0]=="gri":
                constp.update({'rho': node['rho']})
            else:
                constp.update({'rho': params['rho'][3].astype(np.float64)})
        else:
            fitp.update({'rho': params['rho'][3].astype(np.float64)})
            result = np.append(result, fitp['rho'])
        if params['gamma'][0]!="fit":
            if params['gamma'][0]=="gri":
                constp.update({'gamma': node['gamma']})
            else:
                constp.update({'gamma': params['gamma'][3].astype(np.float64)})
        else:
            fitp.update({'gamma': params['gamma'][3].astype(np.float64)})
            result = np.append(result, fitp['gamma'])
        if params['piEE'][0]!="fit":
            if params['piEE'][0]=="gri":
                constp.update({'piEE': node['piEE']})
            else:
                constp.update({'piEE': params['piEE'][3].astype(np.float64)})
        else:
            fitp.update({'piEE': params['piEE'][3].astype(np.float64)})
            result = np.append(result, fitp['piEE'])
        if params['piEN'][0]!="fit":
            if params['piEN'][0]=="gri":
                constp.update({'piEN': node['piEN']})
            else:
                constp.update({'piEN': params['piEN'][3].astype(np.float64)})
        else:
            fitp.update({'piEN': params['piEN'][3].astype(np.float64)})
            result = np.append(result, fitp['piEN'])
        if params['s'][0]!="fit":
            if params['s'][0]=="gri":
                constp.update({'s': node['s']})
            else:
                constp.update({'s': params['s'][3].astype(np.float64)})
        else:
            fitp.update({'s': params['s'][3].astype(np.float64)})
            result = np.append(result, fitp['s'])
        if params['q'][0]!="fit":
            if params['q'][0]=="gri":
                constp.update({'q': node['q']})
            else:
                constp.update({'q': params['q'][3].astype(np.float64)})
        else:
            fitp.update({'q': params['q'][3].astype(np.float64)})
            result = np.append(result, fitp['q'])
        if params['alpha'][0]!="fit":
            if params['alpha'][0]=="gri":
                constp.update({'alpha': node['alpha']})
            else:
                constp.update({'alpha': params['alpha'][3].astype(np.float64)})
        else:
            fitp.update({'alpha': params['alpha'][3].astype(np.float64)})
            result = np.append(result, fitp['alpha'])
        if params['dadt'][0]!="fit":
            if params['dadt'][0]=="gri":
                constp.update({'dadt': node['dadt']})
            else:
                constp.update({'dadt': params['dadt'][3].astype(np.float64)})
        else:
            fitp.update({'dadt': params['dadt'][3].astype(np.float64)})
            result = np.append(result, fitp['dadt'])
        if params['dsdt'][0]!="fit":
            if params['dsdt'][0]=="gri":
                constp.update({'dsdt': node['dsdt']})
            else:
                constp.update({'dsdt': params['dsdt'][3].astype(np.float64)})
        else:
            fitp.update({'dsdt': params['dsdt'][3].astype(np.float64)})
            result = np.append(result, fitp['dsdt'])

        constp.update({'tb': cfgsetup.getfloat('Modelling', 'tb')})
        constp.update({'tb': cfgsetup.getfloat('Modelling', 'tp')})

        fitp = pd.Series(fitp)
        constp = pd.Series(constp)

        # Load models of magnification
        # ----------------------------
        model_names = np.unique(model2load)
        models_address = dict()
        for i in range(len(model_names)):
            name = 'muLAn.models.{:s}'.format(model_names[i])
            models_address.update({model_names[i]: importlib.import_module(name)})

        # Save objects before optimization
        # --------------------------------
        # This is done because some optimization packages works faster with
        # global variables.

        if cfgsetup.getboolean('Modelling', 'Fit'):

            args = dict()
            args.update({'data': pd.DataFrame().from_dict(time_serie)})
            args.update({'fit_params': fitp})
            args.update({'const_params': constp})
            args.update({'instruments': constp})
            args.update({'const_params': constp})
            #args.update({'model_library': pd.DataFrame().from_dict(models_address)})

            # Save file in HDF5 format
            fname = f'args.h5'
            for key, val in args.items():
                val.to_hdf(fname, key=key)

            optimization_names = np.array([cfgsetup.get('Modelling', 'Method')])
            optimization_address = list()
            for i in range(len(optimization_names)):
                name = 'muLAn.models.{:s}'.format(optimization_names[i])
                optimization_address.append(importlib.import_module(name))


        # ------------------------------------------------------------------
        #   Explore the parameters space
        # ------------------------------------------------------------------
        if cfgsetup.getboolean('Modelling', 'Fit'):

            text = "Start minimization..."
            communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

            text = optimization_address[0].help()
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

            try:
                optimization_address[0].search(cfgsetup=cfgsetup, models=models_address,
                                 model_param=model_params, time_serie=time_serie, \
                                 model2load=model2load, interpol_method=interpol_method)
                if os.path.exists(fname): os.remove(fname)
            except:
                if os.path.exists(fname): os.remove(fname)

    if os.path.exists(fname): os.remove(fname)

    # ----------------------------------------------------------------------
    #   Without Data
    # ----------------------------------------------------------------------
    if cfgsetup.getboolean('Plotting', 'Data') == False:

        text = "\n Preparation of a simulation without data "
        communicate(cfgsetup, 1, text, opts=[printoption.reverse])

        # Data file names
        data2find = np.array(cfgobs.options('ObservatoriesDetails'))

        # Cross-match with data to use
        data2use = np.array(cfgsetup.options('Observatories'))
        data2use = np.array([data2find[i] for i in range(len(data2find)) if np.any(data2use == data2find[i])])

        observatories = data2use

        # Reference for plotting
        filename = cfgsetup.get('FullPaths', 'Event') + 'setup.ini'
        file = open(filename, 'r')
        for line in file:
            a = line.split(':')[0]
            if a != '':
                a = a.strip().lower()
                if len(np.where(data2use == a)[0]) > 0:
                    cfgsetup.set('Observatories', 'Reference', a)
                    break

        obs_properties = dict({'name': [], 'colour': [], 'filter': [], 'loc': [], 'exclude': [], 'key': []})
        text = 'Simulated light curve from:\n'
        for i in range(len(observatories)):
            table = [a.strip() for a in cfgobs.get('ObservatoriesDetails', observatories[i]).split(',')]

            obs_properties['name'].append(table[0])
            obs_properties['colour'].append(table[1])
            obs_properties['filter'].append(table[2])
            obs_properties['loc'].append(table[3])
            obs_properties['key'].append(observatories[i])

            if (len(table[4:]) > 0):
                if (table[4] != ''):
                    obs_properties['exclude'].append(table[4:])
                else:
                    obs_properties['exclude'].append(None)
            else:
                obs_properties['exclude'].append(None)

            if i!=len(observatories)-1:
                text = text + '    ' + table[0] + '\n'
            if i==len(observatories)-1:
                text = text + '    ' + table[0]

        communicate(cfgsetup, 1, text)

        # Parameters of modelling
        model_params = dict()
        interpol_method = dict()
        for i in range(len(obs_properties['loc'])):
            name = 'Models_' + obs_properties['loc'][i]
            table = np.array([a.split('/')[1].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            model_params.update({name: table})

            table = np.array([a.split('/')[0].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            name2 = 'DateRanges_' + obs_properties['loc'][i]
            model_params.update({name2: table})

            for j in range(len(unpack_options(cfgsetup, 'Modelling', name))):
                a = unpack_options(cfgsetup, 'Modelling', name)[j]
                a = a.split('/')
                if len(a) == 3:
                    name3 = obs_properties['loc'][i] + '#' + a[1] + '#' + a[0]
                    t1 = float(a[0].split('-')[0].strip())
                    t2 = float(a[0].split('-')[1].strip())
                    # print(name3, t1, t2)
                    interpol_method.update({name3: [np.linspace(t1, t2, int(a[2])), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8')]})
                    # print(interpol_method)

        # Ephemeris
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

        if len(obs_properties['loc']) > 1:
            name1 = obs_properties['loc'][np.where(np.array(
                [obs == cfgsetup.get('Observatories', 'Reference').lower()
                 for obs in observatories]) == True)[0][0]]
            name1 = glob.glob(path + name1 + '.*')[0]
        else:
            name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]

        # Load data
        time_serie = dict({})
        time_serie_temp = dict({})
        model2load = np.array([])
        for i in range(len(observatories)):
            if i == 0:
                time_serie.update({'id': []})
                time_serie.update({'dates': []})
                time_serie.update({'magnitude': []})
                time_serie.update({'err_magn': []})
                time_serie.update({'seeing': []})
                time_serie.update({'background': []})
                time_serie.update({'DsN': []})
                time_serie.update({'DsE': []})
                time_serie.update({'model': []})
                time_serie.update({'interpol': []})
                time_serie.update({'obs': np.array([observatories[i]], dtype='S100')})
                time_serie.update({'loc': np.array([obs_properties['loc'][i]], dtype='S100')})
            else:
                time_serie_temp.update({'id': []})
                time_serie_temp.update({'dates': []})
                time_serie_temp.update({'magnitude': []})
                time_serie_temp.update({'err_magn': []})
                time_serie_temp.update({'seeing': []})
                time_serie_temp.update({'background': []})
                time_serie_temp.update({'DsN': []})
                time_serie_temp.update({'DsE': []})
                time_serie_temp.update({'model': []})
                time_serie_temp.update({'interpol': []})
                time_serie_temp.update({'obs': np.array([observatories[i]], dtype='S100')})
                time_serie_temp.update({'loc': np.array([obs_properties['loc'][i]], dtype='S100')})

                time_serie = {key: np.append(np.array(time_serie[key]),
                                             np.array(time_serie_temp[key])) for key in time_serie}

            # Models
            name = 'Models_' + obs_properties['loc'][i]
            models_temp = model_params[name]
            name = 'DateRanges_' + obs_properties['loc'][i]
            dates_temp = model_params[name]

            key = np.array([key for key in interpol_method])
            for j in range(len(models_temp)):
                model2load = np.append(model2load, models_temp[j])
                time_serie.update({'model': np.append(time_serie['model'], models_temp[j])})

            # Decide if method interpolation is used or not.
            key_list = np.array([key for key in interpol_method])
            if len(key_list) > 0:

                # text = "\nInterpolation of the amplification asked."
                # communicate(cfgsetup, 1, text, opts=[printoption.bright])

                for i in range(len(key_list)):
                    loc = key_list[i].split('#')[0]

                    tmin = float(key_list[i].split('#')[2].split('-')[0])
                    tmax = float(key_list[i].split('#')[2].split('-')[1])

                    # Ephemeris
                    path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Data')

                    if len(obs_properties['loc']) > 1:
                        name1 = obs_properties['loc'][np.where(np.array(
                            [obs == cfgsetup.get('Observatories', 'Reference').lower()
                             for obs in observatories]) == True)[0][0]]
                        name1 = glob.glob(path + name1 + '.*')[0]
                    else:
                        name1 = glob.glob(path + obs_properties['loc'][0] + '.*')[0]

                    c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), \
                                      dec=cfgsetup.get('EventDescription', 'DEC'), frame='icrs')
                    # print(c_icrs.transform_to('barycentrictrueecliptic'))
                    l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
                    b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

                    name2 = glob.glob(path + loc + '.*')[0]
                    sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
                        ephemeris.Ds(name1, name2, l, b, cfgsetup.getfloat('Modelling', 'tp'), \
                                     cfgsetup)

                    if name1 != name2:
                        DsN = DsNs
                        DsE = DsEs
                    else:
                        DsN = DsNe
                        DsE = DsEe

                    interpol_method[key_list[i]][1] = DsN(interpol_method[key_list[i]][0])
                    interpol_method[key_list[i]][2] = DsE(interpol_method[key_list[i]][0])

        # ----------------------------------------------------------------------
        #   Load all the models that will be used
        # ----------------------------------------------------------------------
        def prepar_importation(model2load):
            modules = model2load
            modules_mini = np.array([cfgsetup.get('Modelling', 'Method')])
            path = cfgsetup.get('FullPaths', 'Code') + 'packages/'
            modulesloading_file = open(path + 'modulesloading.py', 'w')

            text = "# -*-coding:Utf-8 -*\n"
            text = text + "import sys\n"
            text = text + "sys.path.insert(0, '" + cfgsetup.get('FullPaths', 'Code') \
                   + "models')\n"

            for i in range(len(modules)):
                text = text + "import " + modules[i] + " as " + modules[i]
                text = text + "\n"

            for i in range(len(modules_mini)):
                text = text + "import " + modules_mini[i] + " as " + modules_mini[i]
                text = text + "\n"

            text = text + "def main():\n"
            text = text + "\tmodels_files = [" + modules[0] + "]\n"

            if len(modules) > 1:
                for i in range(len(modules) - 1):
                    text = text + "\tmodels_files.append(" + modules[i + 1] + ")\n"

            text = text + "\tminim_files = [" + modules_mini[0] + "]\n"

            if len(modules_mini) > 1:
                for i in range(len(modules_mini) - 1):
                    text = text + "\tminim_files.append(" + modules_mini[i + 1] + ")\n"

            text = text + "\treturn models_files, minim_files\n"

            modulesloading_file.write(text)
            modulesloading_file.close()


        model2load = np.unique(model2load)
        prepar_importation(model2load)
        sys.path.insert(0, cfgsetup.get('FullPaths', 'Code') + 'packages/')
        import modulesloading as load_modules

        models_modules, minim_algo = load_modules.main()

        models_modules = {model2load[i]: models_modules[i] for i in range(len(model2load))}

    # Sort the MCMC samples and create a summary file for analysis
    # ============================================================

    if 'sortno' in options:
        text = "Post-processing MCMC output..."
        communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

        if cfgsetup.getboolean('Optimization', 'UseBinaryFiles'):
            if options['sortno']:
                summary = iotools.FitResults(cfgsetup, format='h5')
            else:
                summary = iotools.FitResults(cfgsetup, format='ascii')
                summary.remove_duplicates(inplace=True)
                summary.save(format='h5')
        else:
            if options['sortno']:
                summary = iotools.FitResults(cfgsetup, format='h5')
            else:
                summary = iotools.FitResults(cfgsetup, format='ascii')
                summary.remove_duplicates(inplace=True)
                summary.save(format='h5')
                n = np.min([len(summary.samples), 1000])
                summary.save(format='ascii', N=n)

        # Display a preview in the terminal
        cols = ['fullid', 'dchi2', 'tE', 'q', 's', 'rho', 'tS']
        colnames = [
                'Model ID',
                '\u0394\u03C7\u00b2', 
                'tE', 
                'q', 
                'sep',
                'rho',
                'tS',
                ]
        format={'dchi2': '{:.2f}'.format,
                'tE': " {:.1f}".format,
                'q': " {:.2e}".format,
                's': " {:.2e}".format,
                'rho': " {:.2e}".format,
                'tS': " {:.3f}".format,
                }

        txt = "Best models preview"
        communicate(cfgsetup, 3, txt, opts=[printoption.level1], prefix=False, newline=True, tab=False)
        txt = summary.samples.head(5).to_string(columns=cols, header=colnames, index=False, formatters=format, max_rows=15, line_width=79)
        txt = '   {:s}'.format(txt.replace('\n', '\n   '))
        communicate(cfgsetup, 3, txt, opts=False, prefix=False, newline=False, tab=False)

        txt = "Worst models preview"
        communicate(cfgsetup, 3, txt, opts=[printoption.level1], prefix=False, newline=True, tab=False)
        txt = summary.samples.tail(5).to_string(columns=cols, header=colnames, index=False, formatters=format, max_rows=15, line_width=79)
        txt = '   {:s}'.format(txt.replace('\n', '\n   '))
        communicate(cfgsetup, 3, txt, opts=False, prefix=False, newline=False, tab=False)

    # Compute specified model information
    # ===================================
    # Store data in a new class
    data = Data(time_serie)

    # Select two best models
    x = summary.samples[summary.samples['fullid']<4]
    fname = "{:s}/{:s}.h5".format(
            cfgsetup.get('RelativePaths', 'Archives'),
            cfgsetup.get('Controls', 'Archive'))
    lenses = LensModel(archive=fname)
    lenses.compute(data=data, models=x, magnification=True, lib=models_modules, parser=cfgsetup)

    data_to_plot_models = pd.DataFrame()
    data_to_plot_models['dates'] = np.linspace(5740, 5820, 5000)
    name = 'Models_' + obs_properties['loc'][0]
    models_temp = model_params[name]
    name = 'DateRanges_' + obs_properties['loc'][0]
    dates_temp = model_params[name]

    data_to_plot_models['model'] = 'PSPL'

    for j in range(len(models_temp)):
        model2load = np.append(model2load, models_temp[j])
        tmin = float((dates_temp[j]).split('-')[0].strip())
        tmax = float((dates_temp[j]).split('-')[1].strip())

        mask = (data_to_plot_models['dates'] > tmin) & (data_to_plot_models['dates'] <= tmax)
        data_to_plot_models.loc[mask, 'model'] = models_temp[j]

    # Calculations from ephemeris for parallax
    DsN, DsE = ephemeris.dsndse(l, b, cfgsetup, obs_properties, observatories)

    if (DsN == None) | (DsE == None):
        data_to_plot_models['DsN'] = DsN 
        data_to_plot_models['DsE'] = DsE
    else:
        data_to_plot_models['DsN'] = DsN(data_to_plot_models['dates'].values)
        data_to_plot_models['DsE'] = DsE(data_to_plot_models['dates'].values)

    lenses.magnification_model(epochs=data_to_plot_models, models=x, lib=models_modules)

    sys.exit()
    # Next step is to select only mmodels required for this computation. And then add caustic, residuals etc.

    # ----------------------------------------------------------------------
    #   Plots
    # ----------------------------------------------------------------------
    if cfgsetup.get('Plotting', 'flag'):

        text = "Create the plots..."
        communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

        # List of plots
        plots2load_brut = unpack_options(cfgsetup, 'Plotting', 'Type')

        # Load plot types 
        # ---------------
        plots2load = np.unique(plots2load_brut)
        plottypes_list = dict()
        if(len(plots2load)) > 0:
            for i in range(len(plots2load)):
                text = 'muLAn.plottypes.{:s}'.format(plots2load[i])
                importlib.import_module(text)
                plottypes_list.update({plots2load[i]: getattr(mulanplots, plots2load[i])})

        # Recursively run the plots routines
        for i in range(len(plots2load_brut)):

            text = plots2load_brut[i] + " - " + plottypes_list[plots2load_brut[i]].help()
            communicate(cfgsetup, 1, text, opts=False, prefix=False, newline=True, tab=True)

            options = unpack_options(cfgsetup, 'Plotting', 'Options')[i]
            plottypes_list[plots2load_brut[i]].plot(
                    cfgsetup=cfgsetup, models=models_modules,
                    model_param=model_params, time_serie=time_serie,
                    obs_properties=obs_properties, options=options,
                    interpol_method=interpol_method)

if (__name__ == "__main__"):
    pass

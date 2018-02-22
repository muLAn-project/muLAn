# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# Call the modules asked by the user
# ----------------------------------------------------------------------
# Packages
# ----------------------------------------------------------------------
import os
import sys
import glob
import importlib
import numpy as np
import ConfigParser as cp
from astropy.coordinates import SkyCoord
# ----------------------------------------------------------------------
# Local packages
# ----------------------------------------------------------------------
from general_tools import *
import muLAn.models.ephemeris as ephemeris
import muLAn.models as mulanmodels
import muLAn.plottypes as mulanplots
import muLAn.packages.sortmodels as mulansort

# ====================================================================
# Fonctions
# ====================================================================
def run_sequence(path_event, options):

    # Load configuration files
    cfgsetup = cp.SafeConfigParser()
    cfgsetup.read([path_event + 'setup.ini', path_event + 'advancedsetup.ini'])

    text = "Load parameter files..."
    communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True)

    cfgobs = cp.SafeConfigParser()
    cfgobs.read(path_event + 'observatories.ini')

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

    # Take into account manual options
    if options['plot'] != None:
        cfgsetup.set('Controls', 'Modes', 'Plot')
        cfgsetup.set('Plotting', 'Models', options['plot'])
    if options['fit'] != None:
        cfgsetup.set('Controls', 'Modes', 'Fit')
    cond = (options['fit']) and (options['plot'] != None)
    if cond:
        cfgsetup.set('Controls', 'Modes', 'Fit, Plot')
        cfgsetup.set('Plotting', 'Models', options['plot'])
    if options['archive'] != None:
        print options['archive']
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
    # communicate(cfgsetup, 4, text, opts=False, prefix=False, newline=False)

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
    for i in xrange(len(models_temp)):
        text = "Models for {:s} observations:".format(models_temp[i].split("_")[1].upper())
        communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)
        for j in xrange(len(models_temp2[i])):
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
    [communicate(cfgsetup, 3, text + plots_temp[i] + ": " + options_temp[i], opts=False, prefix=False, newline=False, tab=True) for i in xrange(len(plots_temp))]


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
        data2use = np.array([data2find[i] for i in xrange(len(data2find)) if np.any(data2use == data2find[i])])

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
        data_filenames = [data_filenames[i] for i in xrange(len(data_filenames)) if
                          np.any(data2use == observatories[i].rpartition('.')[0].lower())]
        observatories = [ob.rpartition('.')[0].lower() for ob in observatories if
                         np.any(data2use == ob.rpartition('.')[0].lower())]

        obs_properties = dict({'name': [], 'colour': [], 'filter': [], 'gamma': [],
                'loc': [], 'fluxoumag': [], 'exclude': [], 'key': []})
        # text = 'Data from:\n'
        for i in xrange(len(observatories)):
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
        for i in xrange(len(obs_properties['loc'])):
            name = 'Models_' + obs_properties['loc'][i]
            table = np.array([a.split('/')[1].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            model_params.update({name: table})

            table = np.array([a.split('/')[0].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            name2 = 'DateRanges_' + obs_properties['loc'][i]
            model_params.update({name2: table})

            for j in xrange(len(unpack_options(cfgsetup, 'Modelling', name))):
                a = unpack_options(cfgsetup, 'Modelling', name)[j]
                a = a.split('/')
                if len(a) == 3:
                    name3 = obs_properties['loc'][i] + '#' + a[1] + '#' + a[0]
                    t1 = float(a[0].split('-')[0].strip())
                    t2 = float(a[0].split('-')[1].strip())
                    # print name3, t1, t2
                    interpol_method.update({name3: [np.linspace(t1, t2, int(a[2])), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8')]})
                    # print interpol_method

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

        time_serie = dict({})
        time_serie_temp = dict({})
        model2load = np.array([])
        text_nb_data = ""
        for i in xrange(len(observatories)):
            file = np.loadtxt(data_filenames[i], dtype=format, usecols=(0, 1, 2, 3, 4, 5), \
                              skiprows=0, unpack=False)

            # Rescale the error-bars [see Wyrzykowski et al. (2009)]
            gamma = float(unpack_options(cfgsetup, 'Observatories', observatories[i])[0][1:])
            epsilon = float(unpack_options(cfgsetup, 'Observatories', observatories[i])[1][:-1])

            # Ephemeris
            c_icrs = SkyCoord(ra=cfgsetup.get('EventDescription', 'RA'), \
                              dec=cfgsetup.get('EventDescription', 'DEC'), frame='icrs')
            # print c_icrs.transform_to('barycentrictrueecliptic')
            l = c_icrs.transform_to('barycentrictrueecliptic').lon.degree
            b = c_icrs.transform_to('barycentrictrueecliptic').lat.degree

            name2 = glob.glob(path + obs_properties['loc'][i] + '.*')[0]
            sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs = \
                ephemeris.Ds(name1, name2, l, b, cfgsetup.getfloat('Modelling', 'tp'), \
                             cfgsetup)

            if name1 != name2:
                DsN = DsNs
                DsE = DsEs
            else:
                DsN = DsNe
                DsE = DsEe

            # Timeseries
            dates = file['dates']
            cond = file['dates'] > 2450000
            if cond.sum() > 0:
                dates = file['dates'] - 2450000

            if i == 0:
                time_serie.update({'id': file['id']})
                time_serie.update({'dates': dates})
                time_serie.update({'magnitude': file['magnitude']})
                time_serie.update({'err_magn': np.sqrt(np.power(gamma * file['err_magn'], 2) + epsilon ** 2)})
                time_serie.update({'err_magn_orig': file['err_magn']})
                time_serie.update({'seeing': file['seeing']})
                time_serie.update({'background': file['background']})
                time_serie.update({'DsN': np.full(dates.shape[0], 0, dtype='f8')})
                time_serie.update({'DsE': np.full(dates.shape[0], 0, dtype='f8')})
                time_serie.update({'model': np.full(dates.shape[0], '0', dtype='S100')})
                time_serie.update({'interpol': np.full(dates.shape[0], '0', dtype='S100')})
                time_serie.update({'obs': np.full(dates.shape[0], observatories[i], dtype='S100')})
                time_serie.update({'gamma': np.full(dates.shape[0], obs_properties['gamma'][i], dtype='f8')})
                time_serie.update({'loc': np.full(dates.shape[0], obs_properties['loc'][i], dtype='S100')})
            else:
                time_serie_temp.update({'id': file['id']})
                time_serie_temp.update({'dates': dates})
                time_serie_temp.update({'magnitude': file['magnitude']})
                time_serie_temp.update({'err_magn': np.sqrt(np.power(gamma * file['err_magn'], 2) + epsilon ** 2)})
                time_serie_temp.update({'err_magn_orig': file['err_magn']})
                time_serie_temp.update({'seeing': file['seeing']})
                time_serie_temp.update({'background': file['background']})
                time_serie_temp.update({'DsN': np.full(dates.shape[0], 0, dtype='f8')})
                time_serie_temp.update({'DsE': np.full(dates.shape[0], 0, dtype='f8')})
                time_serie_temp.update({'model': np.full(dates.shape[0], '0', dtype='S100')})
                time_serie_temp.update({'interpol': np.full(dates.shape[0], '0', dtype='S100')})
                time_serie_temp.update({'obs': np.full(dates.shape[0], observatories[i], dtype='S100')})
                time_serie_temp.update({'gamma': np.full(dates.shape[0], obs_properties['gamma'][i], dtype='f8')})
                time_serie_temp.update({'loc': np.full(dates.shape[0], obs_properties['loc'][i], dtype='S100')})

                time_serie = {key: np.append(np.array(time_serie[key]),
                                             np.array(time_serie_temp[key])) for key in time_serie}

            # Select data in the right dates range.
            prop_temp = cfgsetup.get('Observatories', observatories[i]).split(',')

            if len(prop_temp) < 2:
                sys.exit("Error Syntax in [Observatories] from 'setup.ini'.")
            prop_temp = np.delete(prop_temp, [0, 1])
            if len(prop_temp) == 0: prop_temp = ['']

            cond1 = time_serie['obs'] == observatories[i]
            cond_obs = np.invert(np.array(cond1))
            if (len(prop_temp) > 0) & (prop_temp != ['']):
                list_temp = [prop_temp[iii].strip() for iii in xrange(len(prop_temp)) if (prop_temp[iii].strip() != '')]
                for a in list_temp:
                    toinclude = a.split('-')
                    if len(toinclude) == 2:
                        if toinclude[0].strip() == '': toinclude[0] = -999999999.0
                        if toinclude[1].strip() == '': toinclude[1] = 999999999.0
                        cond_inc = np.array(
                            (time_serie['dates'] >= float(toinclude[0])) & (time_serie['dates'] <= float(toinclude[1]))
                            & cond1)
                        cond_obs[cond_inc] = True
                    else:
                        text = 'Problem in the time intervals in setup.ini.'
                        sys.exit(text)

                time_serie = dict({key: time_serie[key][cond_obs] for key in time_serie})

            nb_data_init = (time_serie['obs'] == observatories[i]).sum()

            # Remove data specified in the Observatories files using data IDs.
            prop_temp = cfgobs.get('ObservatoriesDetails', observatories[i]).split(',')
            if len(prop_temp) > 5:
                if prop_temp[5].strip() != '':
                    list_temp = [prop_temp[iii].strip() for iii in xrange(len(prop_temp)) if
                                 (iii > 4) & (prop_temp[iii].strip() != '')]
                    for a in list_temp:
                        toremove = a.split('-')
                        if len(toremove) == 2:
                            cond_rm = np.invert(
                                np.array((time_serie['id'] >= long(toremove[0])) & (time_serie['id'] <= long(toremove[1]))
                                         & (time_serie['obs'] == observatories[i])))
                            time_serie = dict({key: time_serie[key][cond_rm] for key in time_serie})
                        else:
                            cond_rm = np.invert(np.array((time_serie['id'] == long(toremove[0]))
                                                         & (time_serie['obs'] == observatories[i])))
                            time_serie = dict({key: time_serie[key][cond_rm] for key in time_serie})

            a = float(unpack_options(cfgsetup, "Observatories", observatories[i])[2].split("-")[0])
            b = float(unpack_options(cfgsetup, "Observatories", observatories[i])[2].split("-")[1])
            text = "Reading data within {:.6f} --> {:.6f} from {:s}".format(a, b, data_filenames[i].split("/")[-1])
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)
            a = (time_serie['obs'] == observatories[i]).sum()
            b = nb_data_init - a
            if a!=a+b:
                text_nb_data = text_nb_data + "    \033[1m{:6d}\033[0m data + \033[40m\033[97m\033[1m{:6d} data excluded\033[0m\n".format(a, b)
            else:
                text_nb_data = text_nb_data + "    \033[1m{:6d}\033[0m data\n".format(a, b)
            if i==len(observatories)-1:
                text_nb_data = text_nb_data + "  = \033[1m{:6d}\033[0m data in total".format(len(time_serie['dates']))

            # try:
            #     del prop_temp
            #     del list_temp
            #     del toremove
            #     #del cond_rm
            # except:
            #     pass

            # Calculation of Ds
            cond1 = time_serie['obs'] == observatories[i]
            time_serie['DsN'][cond1] = DsN(time_serie['dates'][cond1])
            time_serie['DsE'][cond1] = DsE(time_serie['dates'][cond1])

            # try:
            #     del prop_temp
            #     del cond1
            #     #del cond_obs
            #     del list_temp
            #     del toinclude
            #     del cond_inc
            # except:
            #     pass

            # Models
            name = 'Models_' + obs_properties['loc'][i]
            models_temp = model_params[name]
            name = 'DateRanges_' + obs_properties['loc'][i]
            dates_temp = model_params[name]

            # print time_serie['dates'].shape, time_serie['magnitude'].shape, time_serie['obs'].shape

            key = np.array([key for key in interpol_method])
            for j in xrange(len(models_temp)):
                model2load = np.append(model2load, models_temp[j])
                tmin = float((dates_temp[j]).split('-')[0].strip())
                tmax = float((dates_temp[j]).split('-')[1].strip())

                cond = (time_serie['obs'] == observatories[i]) \
                       & (time_serie['dates'] <= tmax) & (time_serie['dates'] >= tmin)

                time_serie['model'][cond] = models_temp[j]

            cond = time_serie['model'] == '0'
            if cond.sum() > 0:
                time_serie['model'][cond] = models_temp[0]

        communicate(cfgsetup, 3, text_nb_data, opts=False, prefix=False, newline=True, tab=False)

        # del file, tmin, tmax, cond, models_temp, dates_temp, dates
        # del sTe, sEe, sNe, DsTe, DsEe, DsNe, sTs, sEs, sNs, DsTs, DsEs, DsNs
        # del name1, name2, name, format, time_serie_temp

        time_serie.update({'flux': np.power(10, 0.4*(18.0-time_serie['magnitude']))})
        time_serie.update({'err_flux': np.abs((np.log(10) / 2.5) * time_serie['err_magn'] * time_serie['flux'])})
        time_serie.update({'err_flux_orig': np.abs((np.log(10) / 2.5) * time_serie['err_magn_orig'] * time_serie['flux'])})
        time_serie.update({'amp': np.full(time_serie['dates'].shape[0], -1, dtype='f8')})
        time_serie.update({'fs': np.full(time_serie['dates'].shape[0], -999, dtype='f8')})
        time_serie.update({'fb': np.full(time_serie['dates'].shape[0], -999, dtype='f8')})

        # Load flux if available
        # ----------------------
        list_flux = [a for a in observatories
                       for i in xrange(len(obs_properties['key']))
                       if ((a == obs_properties['key'][i])
                           & ((obs_properties['fluxoumag'][i]).lower() == "flux"))]
        for a in list_flux:
            cond = np.where(time_serie['obs'] == a)
            fname = '{:s}{:s}*.dat'.format(
                cfgsetup.get("FullPaths", "Event"),
                cfgsetup.get("RelativePaths", "Data"))
            fname = [aa for aa in glob.glob(fname)
                    if os.path.basename(aa).lower() == '{:s}.dat'.format(a)]
            if len(fname) == 1:
                fname = fname[0]
            else:
                text = "File error when looking for flux."
                sys.exit(text)

            format = {'names': ('id', 'date', 'flux', 'err_flux', 'seeing', 'background'),
                      'formats': ('i8', 'f8', 'f8', 'f8', 'f8', 'f8')}
            file = np.loadtxt(fname, dtype=format, usecols=(0, 1, 2, 3, 4, 5), skiprows=0, unpack=False)
            cond2 = file['date'] > 2450000
            for i in cond[0]:
                if cond2.sum() > 0:
                    idx = (np.abs(file['date'] - time_serie['dates'][i] - 2450000.0)).argmin()
                else:
                    idx = (np.abs(file['date'] - time_serie['dates'][i])).argmin()
                time_serie['flux'][i] = file['flux'][idx]
                time_serie['err_flux_orig'][i] = file['err_flux'][idx]

                # Rescale the error-bars [see Wyrzykowski et al. (2009)]
                gamma = float(unpack_options(cfgsetup, 'Observatories', a)[0][1:])
                epsilon = float(unpack_options(cfgsetup, 'Observatories', a)[1][:-1])

                time_serie['err_flux'][i] = np.sqrt(np.power(gamma * file['err_flux'][idx], 2) + np.power((np.log(10) * epsilon * file['flux'][idx]) / 2.5, 2))

                try:
                    time_serie['magnitude'][i] = 18.0 - 2.5 * np.log10(file['flux'][idx])
                    time_serie['err_magn_orig'][i] = np.abs(2.5 * file['err_flux'][idx] / (file['flux'][idx] * np.log(10)))
                    time_serie['err_magn'][i] = np.sqrt(np.power(gamma * time_serie['err_magn_orig'][i], 2) + epsilon ** 2)
                except:
                    time_serie['magnitude'][i] = 0.0
                    time_serie['err_magn_orig'][i] = 0.0
                    time_serie['err_magn'][i] = 0.0

        # print time_serie['dates'].shape, time_serie['magnitude'].shape, time_serie['obs'].shape

        # Decide if method interpolation is used or not.
        key_list = np.array([key for key in interpol_method])
        if len(key_list) > 0:

            text = "Ask interpolation"
            communicate(cfgsetup, 3, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

            for i in xrange(len(key_list)):
                loc = key_list[i].split('#')[0]

                tmin = float(key_list[i].split('#')[2].split('-')[0])
                tmax = float(key_list[i].split('#')[2].split('-')[1])

                cond1 = (time_serie['dates'] <= tmax) & (time_serie['dates'] >= tmin) &\
                        (time_serie['loc'] == loc)

                # print loc, len(time_serie['model'][cond1]), len(interpol_method[key_list[i]][0])

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
                    # print c_icrs.transform_to('barycentrictrueecliptic')
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

        # Load models of magnification
        # ----------------------------
        models_modules = dict()
        if(len(model2load)) > 0:
            for i in xrange(len(model2load)):
                text = 'muLAn.models.{:s}'.format(model2load[i])
                importlib.import_module(text)
                models_modules.update({model2load[i]: getattr(mulanmodels, model2load[i])})

        # Load minimization algorithms
        # ----------------------------
        modules_mini = np.array([cfgsetup.get('Modelling', 'Method')])
        minim_algo = np.array([])
        if len(modules_mini) > 0:
            for i in xrange(len(modules_mini)):
                text = 'muLAn.models.{:s}'.format(modules_mini[i])
                importlib.import_module(text)
                minim_algo = np.append(minim_algo, getattr(mulanmodels, modules_mini[i]))

        # ------------------------------------------------------------------
        #   Explore the parameters space
        # ------------------------------------------------------------------
        if cfgsetup.getboolean('Modelling', 'Fit'):

            text = "Start minimization..."
            communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)

            text = minim_algo[0].help()
            communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=True, tab=True)

            minim_algo[0].search(cfgsetup=cfgsetup, models=models_modules,
                                 model_param=model_params, time_serie=time_serie, \
                                 model2load=model2load, interpol_method=interpol_method)

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
        data2use = np.array([data2find[i] for i in xrange(len(data2find)) if np.any(data2use == data2find[i])])

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
        for i in xrange(len(observatories)):
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
        for i in xrange(len(obs_properties['loc'])):
            name = 'Models_' + obs_properties['loc'][i]
            table = np.array([a.split('/')[1].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            model_params.update({name: table})

            table = np.array([a.split('/')[0].strip() for a in unpack_options(cfgsetup, 'Modelling', name)])
            name2 = 'DateRanges_' + obs_properties['loc'][i]
            model_params.update({name2: table})

            for j in xrange(len(unpack_options(cfgsetup, 'Modelling', name))):
                a = unpack_options(cfgsetup, 'Modelling', name)[j]
                a = a.split('/')
                if len(a) == 3:
                    name3 = obs_properties['loc'][i] + '#' + a[1] + '#' + a[0]
                    t1 = float(a[0].split('-')[0].strip())
                    t2 = float(a[0].split('-')[1].strip())
                    # print name3, t1, t2
                    interpol_method.update({name3: [np.linspace(t1, t2, int(a[2])), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8'), np.full(int(a[2]), 0, dtype='f8')]})
                    # print interpol_method

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
        for i in xrange(len(observatories)):
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
            for j in xrange(len(models_temp)):
                model2load = np.append(model2load, models_temp[j])
                time_serie.update({'model': np.append(time_serie['model'], models_temp[j])})

            # Decide if method interpolation is used or not.
            key_list = np.array([key for key in interpol_method])
            if len(key_list) > 0:

                # text = "\nInterpolation of the amplification asked."
                # communicate(cfgsetup, 1, text, opts=[printoption.bright])

                for i in xrange(len(key_list)):
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
                    # print c_icrs.transform_to('barycentrictrueecliptic')
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

            for i in xrange(len(modules)):
                text = text + "import " + modules[i] + " as " + modules[i]
                text = text + "\n"

            for i in xrange(len(modules_mini)):
                text = text + "import " + modules_mini[i] + " as " + modules_mini[i]
                text = text + "\n"

            text = text + "def main():\n"
            text = text + "\tmodels_files = [" + modules[0] + "]\n"

            if len(modules) > 1:
                for i in xrange(len(modules) - 1):
                    text = text + "\tmodels_files.append(" + modules[i + 1] + ")\n"

            text = text + "\tminim_files = [" + modules_mini[0] + "]\n"

            if len(modules_mini) > 1:
                for i in xrange(len(modules_mini) - 1):
                    text = text + "\tminim_files.append(" + modules_mini[i + 1] + ")\n"

            text = text + "\treturn models_files, minim_files\n"

            modulesloading_file.write(text)
            modulesloading_file.close()


        model2load = np.unique(model2load)
        prepar_importation(model2load)
        sys.path.insert(0, cfgsetup.get('FullPaths', 'Code') + 'packages/')
        import modulesloading as load_modules

        models_modules, minim_algo = load_modules.main()

        models_modules = {model2load[i]: models_modules[i] for i in xrange(len(model2load))}

    # ----------------------------------------------------------------------
    #   Order models
    # ----------------------------------------------------------------------
    if 'sortno' in options:
        if not options['sortno']:
            text = "Post-process the output files..."
            communicate(cfgsetup, 1, text, opts=[printoption.level0], prefix=True, newline=True, tab=False)
            mulansort.McmcFiles().sort(cfgsetup=cfgsetup)

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
            for i in xrange(len(plots2load)):
                text = 'muLAn.plottypes.{:s}'.format(plots2load[i])
                importlib.import_module(text)
                plottypes_list.update({plots2load[i]: getattr(mulanplots, plots2load[i])})

        # Recursively run the plots routines
        for i in xrange(len(plots2load_brut)):

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

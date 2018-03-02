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

#filename = full_path + '../' + '.pythonexternallibpath'
#file = open(filename, 'r')
#for line in file:
#    path_lib_ext=line
#file.close()
#if path_lib_ext != 'None':
#    sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#    Packages
# ----------------------------------------------------------------------
import os
import sys
import copy
import emcee
import pickle
import glob
import shutil
import datetime
import importlib
import subprocess
import numpy as np
from scipy import stats
from scipy import interpolate
from sklearn import linear_model
import muLAn.models as mulanmodels
# ----------------------------------------------------------------------
#   CLASS
# ----------------------------------------------------------------------
class printoption:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'

    reset = '\033[0m'
    bright = '\033[1m'
    dim = '\033[2m'
    underscore = '\033[4m'
    blink = '\033[5m'
    reverse = '\033[7m'
    hidden = '\033[8m'

    level0 = "\033[1m\033[31m"
    level1 = "\033[1m"
    good = "\033[32m"
# ----------------------------------------------------------------------
#    Functions
# ----------------------------------------------------------------------
def help():
    text = "grid_dmcmc - Differential Markov Chains Monte Carlo."
    return text
# ----------------------------------------------------------------------
def bash_command(text):
    proc = subprocess.Popen(text, shell=True, executable="/bin/bash")
    proc.wait()
# ----------------------------------------------------------------------
# def update_progress(job_title, progress):
def update_progress(job_title, a, b):
    length = 20
    progress = float(a)/float(b)
    block = int(round(length*progress))
    msg = "\r    {0}: [{1}] {2:3.0f}% --> {3:d} / {4:d}".format(job_title, "*"*block + "-"*(length-block), round(progress*100, 2), a, b)
    if progress >= 1: msg = msg + " \033[1m\033[32mDONE\033[0m\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()
# ----------------------------------------------------------------------
def update_progress_grid(a, b, c, d):
    length = 10
    progress = float(a)/float(b)
    progress_grid = float(c) / float(d)
    block = int(round(length*progress))
    block_grid = int(round(length * progress_grid))
    msg = "\r    Grid\033[34m+MCMC\033[0m: {2:4d} / {3:4d} <-- {4:3.0f}% [{0}]\033[34m[{1}] {5:3.0f}% --> {6:d} / {7:d}\033[0m\r".format("*"*block_grid + "-"*(length-block_grid),\
            "*" * block + "-" * (length - block), c, d, round(progress_grid * 100, 2), round(progress*100, 2), a, b)
    sys.stdout.write(msg)
    sys.stdout.flush()
# ----------------------------------------------------------------------
def communicate(cfg, verbose, text, opts=False, prefix=False, newline=False, tab=False):
    if cfg.getint('Modelling', 'Verbose') >= verbose:
        if prefix:
            text = "[muLAn] " + text
        if opts!=False:
            text2=''
            for a in opts:
                text2 = text2 + a
                text = text2 + text + printoption.reset
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print text
        else:
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print text
# ----------------------------------------------------------------------
def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file))
# ----------------------------------------------------------------------
def combin(p):
    Np = len(p)
    Ng = 1
    for m in range(0,Np):
        Ng *= len(p[m])
    gridlist = np.zeros((Np,Ng), dtype='f8') # to test
#    gridlist = np.zeros((Np,Ng))
    Nr = Ng
    for m in range(0,Np):
        Nr = Nr/len(p[m])
        q = 0
        l = 0
        for k in range(0,Ng/Nr):
            for n in range(0,Nr):
                gridlist[m][q] = p[m][l%len(p[m])]
                q += 1
            l += 1
    return gridlist
# ----------------------------------------------------------------------
def binrot(alpha, tau, beta, s, q):
    """Source position over the time. The general conventions used are
    the same as in appendix A in [1].

    Arguments:
    alpha -- the angle between the lens symmetry axis and the source
             trajectory;
    tau -- the time component of the source position;
    beta -- the component of the source position orthogonal to tau;
    s -- the primary-secondary distance;
    q -- the secondary-primary mass ratio.

    Returns:
    x -- numpy array including the x component of the source relative to
         the lens center of mass (CM);
    y -- numpy array including the y component of the source relative to
         the CM.

    References:
    [1] Skowron et al. 2011, 738, 87.
    """

    tau_chap = np.array([np.cos(alpha), np.sin(alpha)])
    beta_chap = np.array([-np.sin(alpha), np.cos(alpha)])

    lenssource = np.array([tau[i] * tau_chap + beta[i] * beta_chap for i in xrange(len(tau))])

    gl1 = s * q/(1+q) * np.array([1, 0])
    lenssource = lenssource - gl1

    return lenssource.T[0], lenssource.T[1]
# ----------------------------------------------------------------------
def test_blending(mb_lim, g_lim, fs, fb, time_serie, cond2):
    g_mod = fb/fs
    mb_mod = 18.0 - 2.5*np.log10(fs+fb)

    # Blending ok
    cond = (g_lim[0] < g_mod) & (g_mod < g_lim[1]) & (mb_lim[0] < mb_mod) & (mb_mod < mb_lim[1])
    if cond:
        fs_new = fs
        fb_new = fb

    # Corners C1
    cond = (g_lim[1] < g_mod) & (mb_lim[1] < mb_mod)
    if cond:
        fs_new = (10**((18.0-mb_lim[1])/2.5))/(1+g_lim[1])
        fb_new = g_lim[1] * fs_new

    # Corners C2
    cond = (g_lim[1] < g_mod) & (mb_mod < mb_lim[0])
    if cond:
        fs_new = (10**(18.0-mb_lim[0]))/(1+g_lim[1])
        fb_new = g_lim[1] * fs_new

    # Corners C3
    cond = (g_mod < g_lim[0]) & (mb_mod < mb_lim[0])
    if cond:
        fs_new = (10**(18.0-mb_lim[0]))/(1+g_lim[0])
        fb_new = g_lim[0] * fs_new

    # Corners C4
    cond = (g_mod < g_lim[0]) & (mb_lim[1] < mb_mod)
    if cond:
        fs_new = (10**((18.0-mb_lim[1])/2.5))/(1+g_lim[0])
        fb_new = g_lim[0] * fs_new

    # Boundary B1
    cond = (g_lim[0] < g_mod) & (g_mod < g_lim[1]) & (mb_lim[1] < mb_mod)
    if cond:
        x = np.atleast_2d(time_serie['amp'][cond2] - 1.0)
        y = np.atleast_2d(time_serie['flux'][cond2] - 10 ** ((18.0 - mb_lim[1]) / 2.5))

        regr = linear_model.LinearRegression(fit_intercept=False)
        regr.fit(x, y)
        fs_new = regr.coef_[0][0]
        fb_new = 10 ** ((18.0 - mb_lim[1]) / 2.5) - fs_new

    # Boundary B2
    cond = (g_lim[1] < g_mod) & (mb_lim[0] < mb_mod) & (mb_mod < mb_lim[1])
    if cond:
        x = np.atleast_2d(time_serie['amp'][cond2] + g_lim[1]).T
        y = np.atleast_2d(time_serie['flux'][cond2]).T

        regr = linear_model.LinearRegression(fit_intercept=False)
        regr.fit(x, y)
        fs_new = regr.coef_[0][0]
        fb_new = g_lim[1]*fs_new

    # Boundary B3
    cond = (g_lim[0] < g_mod) & (g_mod < g_lim[1]) & (mb_mod < mb_lim[0])
    if cond:
        x = np.atleast_2d(time_serie['amp'][cond2] - 1.0)
        y = np.atleast_2d(time_serie['flux'][cond2] - 10 ** ((18.0 - mb_lim[0]) / 2.5))

        regr = linear_model.LinearRegression(fit_intercept=False)
        regr.fit(x, y)
        fs_new = regr.coef_[0][0]
        fb_new = 10 ** ((18.0 - mb_lim[0]) / 2.5) - fs_new

    # Boundary B4
    cond = (g_mod < g_lim[0]) & (mb_lim[0] < mb_mod) & (mb_mod < mb_lim[1])
    if cond:
        x = np.atleast_2d(time_serie['amp'][cond2] + g_lim[0]).T
        y = np.atleast_2d(time_serie['flux'][cond2]).T

        regr = linear_model.LinearRegression(fit_intercept=False)
        regr.fit(x, y)
        fs_new = regr.coef_[0][0]
        fb_new = g_lim[0]*fs_new

    return fs_new, fb_new

# ----------------------------------------------------------------------
# def sort_on_runtime(pos):
#     print pos
#     p = np.atleast_2d(pos)
#     idx = np.argsort(p[:, 3])#[::-1]
#     #print idx
#     return p[idx], idx
# ----------------------------------------------------------------------
def lnprior(param_model):
    p = 0
    if param_model['t0'] < 0:
        p = -np.inf
    if param_model['rho'] < 0:
        p = -np.inf
    if param_model['rho'] > 1.0:
        p = -np.inf
    if param_model['tE'] < 1e-10:
        p = -np.inf
    if param_model['q'] < 1e-9:
        p = -np.inf
    if param_model['q'] > 1.0:
        p = -np.inf
    if param_model['s'] < 1e-10:
        p = -np.inf
    if param_model['s'] > 10:
        p = -np.inf
    return p
# ----------------------------------------------------------------------
def lnprob(theta, time_serie, model_params, fitted_param, nuisance, models_names,
           interpol_method, tb, cfgsetup):
    # print theta[2]

    #import modulesloading as load_modules

    #models, y = load_modules.main()
    #models = {models_names[i] : models[i] for i in xrange(len(models_names))}

    models = dict()
    for i in xrange(len(models_names)):
        text = 'muLAn.models.{:s}'.format(models_names[i])
        importlib.import_module(text)
        models.update({models_names[i]: getattr(mulanmodels, models_names[i])})

    # print res
    # -----------

    # print models
    # print models_names
    # sys.exit()

    # print models
    # print "Hello, c'est moi"

    flag_fix_gamma = 1

    key_list = np.array([])
    for key, value in fitted_param.iteritems():
        key_list = np.append(key_list, key)

    param_model = nuisance
    id=0
    cond = (key_list=='t0')
    if cond.sum()==1:
        param_model.update({'t0' : theta[id]})
        id=id+1
    cond = (key_list=='u0')
    if cond.sum()==1:
        param_model.update({'u0' : theta[id]})
        id=id+1
    cond = (key_list=='tE')
    if cond.sum()==1:
        param_model.update({'tE' : theta[id]})
        id=id+1
    cond = (key_list=='rho')
    if cond.sum()==1:
        param_model.update({'rho' : theta[id]})
        id=id+1
    cond = (key_list=='gamma')
    if cond.sum()==1:
        param_model.update({'gamma' : theta[id]})
        flag_fix_gamma = 0
        id=id+1
    cond = (key_list=='piEE')
    if cond.sum()==1:
        param_model.update({'piEE' : theta[id]})
        id=id+1
    cond = (key_list=='piEN')
    if cond.sum()==1:
        param_model.update({'piEN' : theta[id]})
        id=id+1
    cond = (key_list=='s')
    if cond.sum()==1:
        param_model.update({'s' : theta[id]})
        id=id+1
    cond = (key_list=='q')
    if cond.sum()==1:
        param_model.update({'q' : theta[id]})
        id=id+1
    cond = (key_list=='alpha')
    if cond.sum()==1:
        param_model.update({'alpha' : theta[id]})
        id=id+1
    cond = (key_list == 'dalpha')
    if cond.sum()==1:
        param_model.update({'dalpha' : theta[id]})
        id=id+1
    cond = (key_list == 'ds')
    if cond.sum()==1:
        param_model.update({'ds' : theta[id]})
        # id=id+1

    # Evaluate priors
    chi2 = 0
    lnprior_curr = lnprior(param_model)
    if lnprior_curr != -np.inf:
        # print "Amplification, tu veux ?"
        # Calculation of the amplification
        observatories = np.unique(time_serie['obs'])
        models_lib = np.unique(time_serie['model'])
        for j in xrange(len(observatories)):
            cond2 = (time_serie['obs']==observatories[j])
            #print observatories[j]

            if flag_fix_gamma:
                param_model.update({'gamma': time_serie['gamma'][cond2][0]})

            for i in xrange(models_lib.shape[0]):
                #print models_lib[i]
                cond = (time_serie['model'] == models_lib[i]) & (time_serie['obs']==observatories[j])\
                        & (time_serie['interpol'] == '0')

                if cond.sum() > 0:
                    time_serie_export = time_serie['dates'][cond]
                    DsN_export = time_serie['DsN'][cond]
                    DsE_export = time_serie['DsE'][cond]

                    Ds_export = dict({'N':DsN_export, 'E':DsE_export})

                    try:
                        kwargs_method = dict(cfgsetup.items(models_lib[i]))
                    except:
                        kwargs_method = dict()

                    amp = models[models_lib[i]].magnifcalc(time_serie_export, param_model, Ds=Ds_export, tb=tb, **kwargs_method)

                    time_serie['amp'][cond] = amp

                    del amp

        # Interpolation method
        # -------------------------------------------------------------------------
        key_list = [key for key in interpol_method]

        if len(key_list) > 0:
            for i in xrange(len(key_list)):
                time_serie_export = interpol_method[key_list[i]][0]

                DsN_export = interpol_method[key_list[i]][1]
                DsE_export = interpol_method[key_list[i]][2]

                Ds_export = dict({'N':DsN_export, 'E':DsE_export})

                name = key_list[i].split('#')[1]
                try:
                    kwargs_method = dict(cfgsetup.items(name))
                except:
                    kwargs_method = dict()
                amp = models[name].magnifcalc(time_serie_export, param_model, Ds=Ds_export, tb=tb, **kwargs_method)

                # print amp
                interpol_method[key_list[i]][3] = amp
                interpol_func = interpolate.interp1d(time_serie_export, amp, kind='linear')
                # interpol_func.update({key_list[i]: interpolate.interp1d(time_serie_export, amp)})

                cond = (time_serie['interpol'] == key_list[i])
                if cond.sum() > 0:
                    amp = interpol_func(time_serie['dates'][cond])
                    time_serie['amp'][cond] = amp

        #  Source and blending fluxes.
        # -------------------------------------------------------------------------
        for j in xrange(len(observatories)):
            cond2 = (time_serie['obs']==observatories[j])
            #print observatories[j]

            # Calculation of fs and fb
            # fs, fb = fsfb(time_serie, cond2, blending=True)
            fs, fb = fsfbwsig(time_serie, cond2, blending=True)

            # Relevance of blending for OGLE
            # if (observatories[j]=="ogle-i"):
            #     mb_lim = [17.25, 17.36]
            #     g_lim = [0.0, 10.0]
            #     fs, fb = test_blending(mb_lim, g_lim, fs, fb, time_serie, cond2)

            time_serie['fs'][cond2] = fs
            time_serie['fb'][cond2] = fb

        # print "Amplification, tu as..."
        # Calculation of chi2
        # print param_model, time_serie['amp']

        time_serie['flux_model'] = time_serie['amp']*time_serie['fs'] + time_serie['fb']
        time_serie['chi2pp'] = np.power((time_serie['flux']-time_serie['flux_model'])/time_serie['err_flux'], 2)
        chi2 = np.sum(time_serie['chi2pp'])

    return - chi2/2.0 + lnprior_curr

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
def fsfbwsig(time_serie, cond, blending=True):

    x = np.atleast_2d(time_serie['amp'][cond]).T
    y = np.atleast_2d(time_serie['flux'][cond]).T
    sig = np.atleast_2d(time_serie['err_flux'][cond]).T

    x2 = np.power(x, 2)
    sig2 = np.power(sig, 2)
    s = np.sum(1.0 / sig2)
    sx = np.sum(x / sig2)
    sy = np.sum(y / sig2)
    sxx = np.sum(x2 / sig2)
    sxy = np.sum(x * y / sig2)
    den = s * sxx - sx**2

    if blending:
        fs = (s * sxy - sx * sy) / den
        fb = (sxx * sy - sx * sxy) / den
    else:
        fb = 0.0

    return fs, fb

# ----------------------------------------------------------------------
def ini_chains_gene(fitted_param, nwalkers, params):

    result = []
    key_list = np.array([key for key in fitted_param])
    key_list_order = np.array(['t0', 'u0', 'tE', 'rho', 'gamma', 'piEE', 'piEN', 's', 'q', 'alpha', 'dalpha', 'ds'])
    intersection = np.intersect1d(key_list_order, key_list)
    key_list = [key for key in key_list_order if (len(np.where(intersection==key)[0])>0)]
    l = 0
    while(l<nwalkers):
        table = np.array([])
        for key in key_list:
            if l==0: table = np.append(table, fitted_param[key])
            else:
                a = fitted_param[key] - abs(float(params[key][1]))
                b = fitted_param[key] + abs(float(params[key][2]))
                x = (np.max([a, b]) - np.min([a, b])) * np.random.random_sample() + np.min([a, b])
                table = np.append(table, x)
        result.append(table)
        l = l + 1
    return result

# ----------------------------------------------------------------------
#   Differential MCMC
# ----------------------------------------------------------------------
def search(cfgsetup=False, models=False, model_param=False, time_serie=False,\
           model2load=False, interpol_method=False):

    # ==================================================================
    #   Preparing MCMC
    # ==================================================================
    # Emergency Stop initialization
    if os.path.exists(cfgsetup.get('FullPaths', 'Event') + '.emergencystop'):
        os.remove(cfgsetup.get('FullPaths', 'Event') + '.emergencystop')

    file = open(cfgsetup.get('FullPaths', 'Event') + '.emergencystop', 'w')
    file.write('0')
    file.close()

    fn_lock = cfgsetup.get('FullPaths', 'Event') + '.lock'
    if not os.path.exists(fn_lock): open(fn_lock, 'w').close()

    # Parameter to be fitted / Nuisance parameters
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
        'dalpha': np.array([a.strip() for a in cfgsetup.get('Modelling', 'dalpha').split(',')]),\
        'ds': np.array([a.strip() for a in cfgsetup.get('Modelling', 'ds').split(',')])\
        }
    
    # Files
    path = cfgsetup.get('FullPaths', 'Event')\
            + cfgsetup.get('RelativePaths', 'ModelsHistory')\
            + cfgsetup.get('Controls', 'Archive')\
            + '-ModelsSummary.txt'
    if os.path.exists(path): os.remove(path)

    sys.path.insert(0, cfgsetup.get('FullPaths', 'Code') + 'packages/')
    if (cfgsetup.getboolean('FitSetupDMCMC', 'Resume')==False):

        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')
        shutil.rmtree(path)
        if not os.path.exists(path): os.makedirs(path)

        for i in xrange(cfgsetup.getint('FitSetupDMCMC', 'Chains')):
            filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                + '-c{:04d}'.format(i) + '.txt'
            file_chains = open(filename4chains, 'a')
            format4chains = '# Exploration of chain n°{:d}.\n'.format(i)\
                + '#{:>9}   '.format('ID')\
                + '{:>17}   '.format('t0')\
                + '{:>17}   '.format('u0')\
                + '{:>17}   '.format('tE')\
                + '{:>17}   '.format('rho')\
                + '{:>17}   '.format('gamma')\
                + '{:>17}   '.format('piEN')\
                + '{:>17}   '.format('piEE')\
                + '{:>17}   '.format('s')\
                + '{:>17}   '.format('q')\
                + '{:>17}   '.format('alpha')\
                + '{:>17}   '.format('dalpha')\
                + '{:>17}   '.format('ds')\
                + '{:>17}   '.format('chi2')\
                + '{:>7}   '.format('accrate')\
                + '{:>8}   '.format('date')\
                + '{:>6}   '.format('hour')\
                + '{:>17}   '.format('chi2/dof')\
                + '\n'
            file_chains.write(format4chains)
            file_chains.close()
            accrate_loaded = np.array([])
            id_loaded = np.array([])
    else:
        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')
        if not os.path.exists(path):
            text = "\n\033[1m\033[91mDirectory with chains is missing in 'Resume' mode. muLAn killed.\033[0m"
            sys.exit(text)
        else:
            path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')

            fnames_chains = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*-c*.txt")
            fnames_chains_exclude = glob.glob(path + cfgsetup.get('Controls', 'Archive') + "*g*.txt")

            temp = []
            for a in fnames_chains:
                if (a in fnames_chains_exclude) == False:
                    temp.append(a)

            fnames_chains = copy.deepcopy(temp)
            del temp, fnames_chains_exclude

            nb_chains = len(fnames_chains)

            if nb_chains != cfgsetup.getint("FitSetupDMCMC", "Chains"):
                text = "\n\033[1m\033[91mThe number of chains does not fit in 'Resume' mode. muLAn killed.\033[0m"
                sys.exit(text)

            samples_file = dict({'chi2': [], 't0': [], 'u0': [], 'tE': [], 'rho': [], 'gamma': [], 'piEE': [], 'piEN': [], 's': [], 'q': [], 'alpha': [], 'dalpha': [], 'ds': [], 'chain': [], 'fullid': [], 'date_save': [], 'time_save': [], 'id': [], 'accrate': [], 'chi2/dof': []})

            accrate_loaded = np.array([])
            id_loaded = np.array([])
            for i in xrange(nb_chains):

                file = open(fnames_chains[i], 'r')
                for line in file:
                    params_model = line

                    if params_model[0] == '#':
                        continue

                    samples_file['id'].append(int([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][0]))
                    samples_file['t0'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][1]))
                    samples_file['u0'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][2]))
                    samples_file['tE'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][3]))
                    samples_file['rho'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][4]))
                    samples_file['gamma'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][5]))
                    samples_file['piEN'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][6]))
                    samples_file['piEE'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][7]))
                    samples_file['s'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][8]))
                    samples_file['q'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][9]))
                    samples_file['alpha'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][10]))
                    samples_file['dalpha'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][11]))
                    samples_file['ds'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][12]))
                    samples_file['chi2'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][13]))
                    samples_file['accrate'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][14]))
                    samples_file['date_save'].append(int([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][15]))
                    samples_file['time_save'].append([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][16])
                    samples_file['chi2/dof'].append(float([a for a in (params_model.split('\n')[0].split(' ')) if (a != '')][17]))
                    samples_file['chain'].append(int(fnames_chains[i][-8:-4]))
                    samples_file['fullid'].append(-1)
                file.close()
                accrate_loaded = np.append(accrate_loaded, samples_file['accrate'][-1])
                id_loaded = np.append(id_loaded, samples_file['id'][-1])

            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')\
                  + cfgsetup.get('Controls', 'Archive') + "-lastposition.p"

            file = open(filename, "r")
            pos_pickle = pickle.load(file)
            file.close()

            filename = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')\
                       + cfgsetup.get('Controls', 'Archive') + "-rgenerator.p"

            file = open(filename, "r")
            rgenerator_piclke = pickle.load(file)
            file.close()

            del samples_file

    # Prepare the grids
    grid_params = np.array([])
    format = 'import numpy as np\nimport pickle\ntab = np.array(['
    i = 0
    if params['t0'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['t0'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['t0'][1], params['t0'][2], params['t0'][3])
    if params['u0'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['u0'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['u0'][1], params['u0'][2], params['u0'][3])
    if params['tE'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['tE'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['tE'][1], params['tE'][2], params['tE'][3])
    if params['rho'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['rho'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['rho'][1], params['rho'][2], params['rho'][3])
    if params['gamma'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['gamma'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['gamma'][1], params['gamma'][2], params['gamma'][3])
    if params['piEE'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['piEE'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['piEE'][1], params['piEE'][2], params['piEE'][3])
    if params['piEN'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['piEN'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['piEN'][1], params['piEN'][2], params['piEN'][3])
    if params['s'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['s'])
        a = float(params['s'][1])
        b = float(params['s'][2])
        if (a > 0) & (b > 0):
            format = format + 'np.logspace({0:.10e}, {1:.10e}, {2:s}),'.format(np.log10(a), np.log10(b), params['s'][3])
        else:
            sys.exit('Please enter a positive value for s.')
    if params['q'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['q'])
        a = float(params['q'][1])
        b = float(params['q'][2])
        if (a > 0) & (b > 0):
            format = format + 'np.logspace({0:.10e}, {1:.10e}, {2:s}),'.format(np.log10(a), np.log10(b), params['q'][3])
        else:
            sys.exit('Please enter a positive value for q.')
    if params['alpha'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['alpha'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['alpha'][1], params['alpha'][2], params['alpha'][3])
    if params['dalpha'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['dalpha'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['dalpha'][1], params['dalpha'][2], params['dalpha'][3])
    if params['ds'][0]=="gri":
        i = i + 1
        grid_params = np.append(grid_params, ['ds'])
        format = format + 'np.linspace({0:s}, {1:s}, {2:s}),'.format(params['ds'][1], params['ds'][2], params['ds'][3])

    format = format[:-1] + '])\n'
    format = format + 'file_save = open("' + cfgsetup.get('FullPaths', 'Code')\
        + 'tmp.p", "w")\npickle.dump(tab, file_save)\nfile_save.close()\n'

    filename = cfgsetup.get('FullPaths', 'Code') + 'temp_grid.py'
    file_temp = open(filename, 'w')
    file_temp.write(format)
    file_temp.close()

    flag_grid_yes = 1
    if i>0:
        execfile(filename)

        filename = cfgsetup.get('FullPaths', 'Code') + 'tmp.p'
        file = open(filename, 'r')
        grid_values = pickle.load(file)
        file.close()

        os.remove(filename)

        grid_values_combined = combin(grid_values)

        nb_params_grid = len(grid_values_combined)
        lengh_grid = len(grid_values_combined.T)
    else:
        nb_params_grid = 1
        lengh_grid = 1
        flag_grid_yes = 0

    filename = cfgsetup.get('FullPaths', 'Code') + 'temp_grid.py'
    os.remove(filename)

    # Prepare the DMCMC
    # print lengh_grid
    t0_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    u0_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    tE_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    rho_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    gamma_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    piEN_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    piEE_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    s_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    q_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    alpha_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    dalpha_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    ds_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    lnprob_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    accrate_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'f8')
    date_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'S8')
    hour_best = np.empty(cfgsetup.getint('FitSetupDMCMC', 'Chains'), 'S6')

    for id_grid in xrange(lengh_grid):

        # if flag_grid_yes:
            # text = '\nGrid: {:d} / {:d}'.format(id_grid+1, lengh_grid)
            # communicate(cfgsetup, 1, text)

        # update_progress_grid("Grid progression", id_grid+1, lengh_grid)

        if flag_grid_yes:
            node = dict()
            for id2_grid in xrange(nb_params_grid):
                node.update({grid_params[id2_grid] : grid_values_combined.T[id_grid][id2_grid]})

            path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')

            for i in xrange(cfgsetup.getint('FitSetupDMCMC', 'Chains')):
                filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                    + '-c{:04d}'.format(i) + '-g{:d}'.format(id_grid) + '.txt'
                file_chains = open(filename4chains, 'a')
                format4chains = '# Exploration of chain n°{:d}.\n'.format(i)\
                    + '#{:>9}   '.format('ID')\
                    + '{:>17}   '.format('t0')\
                    + '{:>17}   '.format('u0')\
                    + '{:>17}   '.format('tE')\
                    + '{:>17}   '.format('rho')\
                    + '{:>17}   '.format('gamma')\
                    + '{:>17}   '.format('piEN')\
                    + '{:>17}   '.format('piEE')\
                    + '{:>17}   '.format('s')\
                    + '{:>17}   '.format('q')\
                    + '{:>17}   '.format('alpha')\
                    + '{:>17}   '.format('dalpha')\
                    + '{:>17}   '.format('ds')\
                    + '{:>17}   '.format('chi2')\
                    + '{:>7}   '.format('accrate')\
                    + '{:>8}   '.format('date')\
                    + '{:>6}   '.format('hour')\
                    + '{:>17}   '.format('chi2/dof')\
                    + '\n'
                file_chains.write(format4chains)
                file_chains.close()

        # Prepare the DMCMC
        nuisance = dict()
        fitted_param = dict()
        result = np.array([])
        if params['t0'][0]!="fit":
            if params['t0'][0]=="gri":
                nuisance.update({'t0': node['t0']})
            else:
                nuisance.update({'t0': params['t0'][3].astype(np.float64)})
        else:
            fitted_param.update({'t0': params['t0'][3].astype(np.float64)})
            result = np.append(result, fitted_param['t0'])
        if params['u0'][0]!="fit":
            if params['u0'][0]=="gri":
                nuisance.update({'u0': node['u0']})
            else:
                nuisance.update({'u0': params['u0'][3].astype(np.float64)})
        else:
            fitted_param.update({'u0': params['u0'][3].astype(np.float64)})
            result = np.append(result, fitted_param['u0'])
        if params['tE'][0]!="fit":
            if params['tE'][0]=="gri":
                nuisance.update({'tE': node['tE']})
            else:
                nuisance.update({'tE': params['tE'][3].astype(np.float64)})
        else:
            fitted_param.update({'tE': params['tE'][3].astype(np.float64)})
            result = np.append(result, fitted_param['tE'])
        if params['rho'][0]!="fit":
            if params['rho'][0]=="gri":
                nuisance.update({'rho': node['rho']})
            else:
                nuisance.update({'rho': params['rho'][3].astype(np.float64)})
        else:
            fitted_param.update({'rho': params['rho'][3].astype(np.float64)})
            result = np.append(result, fitted_param['rho'])
        if params['gamma'][0]!="fit":
            if params['gamma'][0]=="gri":
                nuisance.update({'gamma': node['gamma']})
            else:
                nuisance.update({'gamma': params['gamma'][3].astype(np.float64)})
        else:
            fitted_param.update({'gamma': params['gamma'][3].astype(np.float64)})
            result = np.append(result, fitted_param['gamma'])
        if params['piEE'][0]!="fit":
            if params['piEE'][0]=="gri":
                nuisance.update({'piEE': node['piEE']})
            else:
                nuisance.update({'piEE': params['piEE'][3].astype(np.float64)})
        else:
            fitted_param.update({'piEE': params['piEE'][3].astype(np.float64)})
            result = np.append(result, fitted_param['piEE'])
        if params['piEN'][0]!="fit":
            if params['piEN'][0]=="gri":
                nuisance.update({'piEN': node['piEN']})
            else:
                nuisance.update({'piEN': params['piEN'][3].astype(np.float64)})
        else:
            fitted_param.update({'piEN': params['piEN'][3].astype(np.float64)})
            result = np.append(result, fitted_param['piEN'])
        if params['s'][0]!="fit":
            if params['s'][0]=="gri":
                nuisance.update({'s': node['s']})
            else:
                nuisance.update({'s': params['s'][3].astype(np.float64)})
        else:
            fitted_param.update({'s': params['s'][3].astype(np.float64)})
            result = np.append(result, fitted_param['s'])
        if params['q'][0]!="fit":
            if params['q'][0]=="gri":
                nuisance.update({'q': node['q']})
            else:
                nuisance.update({'q': params['q'][3].astype(np.float64)})
        else:
            fitted_param.update({'q': params['q'][3].astype(np.float64)})
            result = np.append(result, fitted_param['q'])
        if params['alpha'][0]!="fit":
            if params['alpha'][0]=="gri":
                nuisance.update({'alpha': node['alpha']})
            else:
                nuisance.update({'alpha': params['alpha'][3].astype(np.float64)})
        else:
            fitted_param.update({'alpha': params['alpha'][3].astype(np.float64)})
            result = np.append(result, fitted_param['alpha'])
        if params['dalpha'][0]!="fit":
            if params['dalpha'][0]=="gri":
                nuisance.update({'dalpha': node['dalpha']})
            else:
                nuisance.update({'dalpha': params['dalpha'][3].astype(np.float64)})
        else:
            fitted_param.update({'dalpha': params['dalpha'][3].astype(np.float64)})
            result = np.append(result, fitted_param['dalpha'])
        if params['ds'][0]!="fit":
            if params['ds'][0]=="gri":
                nuisance.update({'ds': node['ds']})
            else:
                nuisance.update({'ds': params['ds'][3].astype(np.float64)})
        else:
            fitted_param.update({'ds': params['ds'][3].astype(np.float64)})
            result = np.append(result, fitted_param['ds'])

        # Parameters of MCMC
        ndim, nwalkers = len(fitted_param), cfgsetup.getint('FitSetupDMCMC', 'Chains')  # Attention nwalkers nombre pair.
        # pos = [result + 0.1*np.random.randn(ndim) for i in xrange(
        #         nwalkers)]
        if fitted_param!={}:
            # Use delta random in a specified interval
            if cfgsetup.getboolean("FitSetupDMCMC", "Resume"):
                if pos_pickle.shape[1] != len(fitted_param):
                    text = "\n\033[1m\033[91mThe number of fitted parameters does not fit in 'Resume' mode. muLAn killed.\033[0m"
                    sys.exit(text)
                if flag_grid_yes==1:
                    text = "\n\033[1m\033[91m'Resume' mode not compatible with a grid. muLAn killed.\033[0m"
                    sys.exit(text)

                pos = pos_pickle
                rstate = rgenerator_piclke
            else:
                pos = ini_chains_gene(fitted_param, nwalkers, params)
                rstate = None

            # Sampler
            if id_grid > 0:
                del sampler
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                    args=(time_serie, model_param, fitted_param, nuisance, model2load, interpol_method, cfgsetup.getfloat('Modelling', 'tb'), cfgsetup),
                    threads=cfgsetup.getint('FitSetupDMCMC', 'Threads'))
            # ==============================================================
            #   RUN MCMC
            # ==============================================================
            chain_lengh = cfgsetup.getint('FitSetupDMCMC', 'ChainLength')

            id_model = np.ones(nwalkers, dtype=np.int64)
            for pos, lnprobc, rstate in sampler.sample(pos, rstate0=rstate, iterations=chain_lengh, storechain=False):

                # text = 'MCMC: {:d} / {:d}'.format(id_model[0], chain_lengh)
                # communicate(cfgsetup, 1, text)

                if cfgsetup.getint("Modelling", "Verbose") >=3:
                    if flag_grid_yes:
                        update_progress_grid(id_model[0], chain_lengh, id_grid+1, lengh_grid)
                    else:
                        update_progress("Progression", id_model[0], chain_lengh)

                accrate = sampler.acceptance_fraction

                key_list = np.array([])
                for key, value in fitted_param.iteritems():
                    key_list = np.append(key_list, key)

                for i in xrange(nwalkers):
                    param_model = nuisance
                    id = 0
                    cond = (key_list=='t0')
                    if cond.sum()==1:
                        param_model.update({'t0' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='u0')
                    if cond.sum()==1:
                        param_model.update({'u0' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='tE')
                    if cond.sum()==1:
                        param_model.update({'tE' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='rho')
                    if cond.sum()==1:
                        param_model.update({'rho' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='gamma')
                    if cond.sum()==1:
                        param_model.update({'gamma' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='piEE')
                    if cond.sum()==1:
                        param_model.update({'piEE' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='piEN')
                    if cond.sum()==1:
                        param_model.update({'piEN' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='s')
                    if cond.sum()==1:
                        param_model.update({'s' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='q')
                    if cond.sum()==1:
                        param_model.update({'q' : pos[i][id]})
                        id=id+1
                    cond = (key_list=='alpha')
                    if cond.sum()==1:
                        param_model.update({'alpha' : pos[i][id]})
                        id=id+1
                    cond = (key_list == 'dalpha')
                    if cond.sum() == 1:
                        param_model.update({'dalpha': pos[i][id]})
                        id=id+1
                    cond = (key_list == 'ds')
                    if cond.sum() == 1:
                        param_model.update({'ds': pos[i][id]})
                        id=id+1
                        # id=id+1

                    if flag_grid_yes:
                        # Best model
                        if id_model[i]==1:
                            t0_best[i] = param_model['t0']
                            u0_best[i] = param_model['u0']
                            tE_best[i] = param_model['tE']
                            rho_best[i] = param_model['rho']
                            gamma_best[i] = param_model['gamma']
                            piEN_best[i] = param_model['piEN']
                            piEE_best[i] = param_model['piEE']
                            s_best[i] = param_model['s']
                            q_best[i] = param_model['q']
                            alpha_best[i] = param_model['alpha']
                            dalpha_best[i] = param_model['dalpha']
                            ds_best[i] = param_model['ds']
                            lnprob_best[i] = lnprobc[i]
                            accrate_best[i] = accrate[i]
                            date_best[i] = datetime.date.today().strftime("%Y%m%d")
                            hour_best[i] = datetime.datetime.utcnow().strftime("%H%M%S")

                        elif lnprobc[i]>lnprob_best[i]:
                            t0_best[i] = param_model['t0']
                            u0_best[i] = param_model['u0']
                            tE_best[i] = param_model['tE']
                            rho_best[i] = param_model['rho']
                            gamma_best[i] = param_model['gamma']
                            piEN_best[i] = param_model['piEN']
                            piEE_best[i] = param_model['piEE']
                            s_best[i] = param_model['s']
                            q_best[i] = param_model['q']
                            alpha_best[i] = param_model['alpha']
                            dalpha_best[i] = param_model['dalpha']
                            ds_best[i] = param_model['ds']
                            lnprob_best[i] = lnprobc[i]
                            accrate_best[i] = accrate[i]
                            date_best[i] = datetime.date.today().strftime("%Y%m%d")
                            hour_best[i] = datetime.datetime.utcnow().strftime("%H%M%S")

                        # Save Chains
                        filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                            + '-c{:04d}'.format(i) + '-g{:d}'.format(id_grid) + '.txt'
                        file_chains = open(filename4chains, 'a')
                        format4chains = '{:>10d}   '.format(id_model[i])\
                            + '{:+.10e}   '.format(param_model['t0'])\
                            + '{:+.10e}   '.format(param_model['u0'])\
                            + '{:+.10e}   '.format(param_model['tE'])\
                            + '{:+.10e}   '.format(param_model['rho'])\
                            + '{:+.10e}   '.format(param_model['gamma'])\
                            + '{:+.10e}   '.format(param_model['piEN'])\
                            + '{:+.10e}   '.format(param_model['piEE'])\
                            + '{:+.10e}   '.format(param_model['s'])\
                            + '{:+.10e}   '.format(param_model['q'])\
                            + '{:+.10e}   '.format(param_model['alpha'])\
                            + '{:+.10e}   '.format(param_model['dalpha'])\
                            + '{:+.10e}   '.format(param_model['ds'])\
                            + '{:+.10e}   '.format(-2.0*lnprobc[i])\
                            + '{:>7.3f}   '.format(accrate[i])\
                            + '{:8}   '.format(datetime.date.today().strftime("%Y%m%d"))\
                            + '{:6}   '.format(datetime.datetime.utcnow().strftime("%H%M%S"))\
                            + '{:+.10e}'.format(-2.0*lnprobc[i]/(len(time_serie['dates'])-len(fitted_param)-len(grid_params)))\
                            + '\n'
                        file_chains.write(format4chains)
                        file_chains.close()
                    else:

                        if (len(accrate_loaded)>0) & (len(id_loaded)>0):
                            id_model_curr = int(id_model[i] + id_loaded[i])
                            accrate_curr = (1.0 * accrate[i] * id_model[i] + 1.0 * accrate_loaded[i] * id_loaded[i]) / id_model_curr
                        else:
                            id_model_curr = int(id_model[i])
                            accrate_curr = accrate[i]

                        path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')
                        filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                            + '-c{:04d}'.format(i) + '.txt'
                        file_chains = open(filename4chains, 'a')
                        format4chains = '{:>10d}   '.format(id_model_curr)\
                            + '{:+.10e}   '.format(param_model['t0'])\
                            + '{:+.10e}   '.format(param_model['u0'])\
                            + '{:+.10e}   '.format(param_model['tE'])\
                            + '{:+.10e}   '.format(param_model['rho'])\
                            + '{:+.10e}   '.format(param_model['gamma'])\
                            + '{:+.10e}   '.format(param_model['piEN'])\
                            + '{:+.10e}   '.format(param_model['piEE'])\
                            + '{:+.10e}   '.format(param_model['s'])\
                            + '{:+.10e}   '.format(param_model['q'])\
                            + '{:+.10e}   '.format(param_model['alpha'])\
                            + '{:+.10e}   '.format(param_model['dalpha'])\
                            + '{:+.10e}   '.format(param_model['ds'])\
                            + '{:+.10e}   '.format(-2.0*lnprobc[i])\
                            + '{:>7.3f}   '.format(accrate_curr)\
                            + '{:8}   '.format(datetime.date.today().strftime("%Y%m%d"))\
                            + '{:6}   '.format(datetime.datetime.utcnow().strftime("%H%M%S"))\
                            + '{:+.10e}'.format(-2.0*lnprobc[i]/(len(time_serie['dates'])-len(fitted_param)-len(grid_params)))\
                            + '\n'
                        file_chains.write(format4chains)
                        file_chains.close()
                    id_model[i] = id_model[i] + 1

                # Emergency Stop
                file = open(cfgsetup.get('FullPaths', 'Event') + '.emergencystop', 'r')
                stop = 0
                for line in file:
                    if line.strip() == '1':
                        stop=1
                file.close()

                fn_lock = cfgsetup.get('FullPaths', 'Event') + '.lock'
                if not os.path.exists(fn_lock): stop=1

                # Record the last position
                filename4pos = path + cfgsetup.get('Controls', 'Archive') + '-lastposition.p'
                file_save = open(filename4pos, "w")
                pickle.dump(pos, file_save)
                file_save.close()
                del pos

                # Record the state of the pseudo-random generator
                filename = path + cfgsetup.get('Controls', 'Archive') + '-rgenerator.p'
                file_save = open(filename, "w")
                pickle.dump(rstate, file_save)
                file_save.close()

                if stop==1:
                    break

                # print "On continue"

            # Save best model if grid
            if flag_grid_yes:
                for i in xrange(nwalkers):
                    path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')
                    filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                        + '-c{:04d}'.format(i) + '.txt'
                    file_chains = open(filename4chains, 'a')
                    format4chains = '{:>10d}   '.format(id_grid+1)\
                        + '{:+.10e}   '.format(t0_best[i])\
                        + '{:+.10e}   '.format(u0_best[i])\
                        + '{:+.10e}   '.format(tE_best[i])\
                        + '{:+.10e}   '.format(rho_best[i])\
                        + '{:+.10e}   '.format(gamma_best[i])\
                        + '{:+.10e}   '.format(piEN_best[i])\
                        + '{:+.10e}   '.format(piEE_best[i])\
                        + '{:+.10e}   '.format(s_best[i])\
                        + '{:+.10e}   '.format(q_best[i])\
                        + '{:+.10e}   '.format(alpha_best[i])\
                        + '{:+.10e}   '.format(dalpha_best[i])\
                        + '{:+.10e}   '.format(ds_best[i]) \
                        + '{:+.10e}   '.format(-2.0*lnprob_best[i])\
                        + '{:>7.3f}   '.format(accrate[i])\
                        + '{:8}   '.format(date_best[i])\
                        + '{:6}   '.format(hour_best[i])\
                        + '{:+.10e}'.format(-2.0*lnprobc[i]/(len(time_serie['dates'])-len(fitted_param)-len(grid_params)))\
                        + '\n'
                    file_chains.write(format4chains)
                    file_chains.close()

            if stop==1:
                break

            # Create an archive for each MCMC on the grid
            if flag_grid_yes:
                path_event = cfgsetup.get('FullPaths', 'Event')

                path = path_event + cfgsetup.get('RelativePaths', 'Chains')\
                    + '{:d}'.format(id_grid) + '/'
                os.makedirs(path)

                text = 'cp ' + path_event + cfgsetup.get('RelativePaths', 'Chains')\
                    + '*g' + '{:d}'.format(id_grid) + '* ' + path
                bash_command(text)

                text = 'cp ' + path_event + cfgsetup.get('RelativePaths', 'Chains')\
                    + '*.p ' + path
                bash_command(text)

                shutil.make_archive(path, 'zip', path)
                shutil.rmtree(path)

                text = 'rm ' + path_event + cfgsetup.get('RelativePaths', 'Chains')\
                    + '*g' + '{:d}'.format(id_grid) + '* '
                bash_command(text)
        else:
            stop = 0
            for i in xrange(nwalkers):
                param_model = nuisance


                # Calculation of the amplification
                observatories = np.unique(time_serie['obs'])
                models_lib = np.unique(time_serie['model'])
                for jjj in xrange(len(observatories)):
                    cond2 = (time_serie['obs']==observatories[jjj])

                    for iii in xrange(models_lib.shape[0]):
                        cond = (time_serie['model'] == models_lib[iii]) & (time_serie['obs']==observatories[jjj])

                        if cond.sum() > 0:
                            time_serie_export = time_serie['dates'][cond]
                            DsN_export = time_serie['DsN'][cond]
                            DsE_export = time_serie['DsE'][cond]

                            Ds_export = dict({'N':DsN_export, 'E':DsE_export})

                            try:
                                kwargs_method = dict(cfgsetup.items(models_lib[iii]))
                            except:
                                kwargs_method = dict()

                            amp = models[models_lib[iii]].magnifcalc(time_serie_export, param_model, Ds=Ds_export, tb=cfgsetup.getfloat('Modelling', 'tb'), **kwargs_method)

                            time_serie['amp'][cond] = amp

                            del amp

                    # Calculation of fs and fb
                    # fs, fb = fsfb(time_serie, cond2, blending=True)
                    fs, fb = fsfbwsig(time_serie, cond2, blending=True)

                    time_serie['fs'][cond2] = fs
                    time_serie['fb'][cond2] = fb

                # Calculation of chi2
                time_serie['flux_model'] = time_serie['amp']*time_serie['fs'] + time_serie['fb']
                time_serie['chi2pp'] = np.power((time_serie['flux']-time_serie['flux_model'])/time_serie['err_flux'], 2)
                chi2_ini = np.sum(time_serie['chi2pp'])

                path = cfgsetup.get('FullPaths', 'Event') + cfgsetup.get('RelativePaths', 'Chains')
                filename4chains = path + cfgsetup.get('Controls', 'Archive')\
                    + '-c{:04d}'.format(i) + '.txt'
                file_chains = open(filename4chains, 'a')
                format4chains = '{:>10d}   '.format(0)\
                    + '{:+.10e}   '.format(param_model['t0'])\
                    + '{:+.10e}   '.format(param_model['u0'])\
                    + '{:+.10e}   '.format(param_model['tE'])\
                    + '{:+.10e}   '.format(param_model['rho'])\
                    + '{:+.10e}   '.format(param_model['gamma'])\
                    + '{:+.10e}   '.format(param_model['piEN'])\
                    + '{:+.10e}   '.format(param_model['piEE'])\
                    + '{:+.10e}   '.format(param_model['s'])\
                    + '{:+.10e}   '.format(param_model['q'])\
                    + '{:+.10e}   '.format(param_model['alpha'])\
                    + '{:+.10e}   '.format(param_model['dalpha'])\
                    + '{:+.10e}   '.format(param_model['ds'])\
                    + '{:+.10e}   '.format(chi2_ini)\
                    + '{:>7.3f}   '.format(0.0)\
                    + '{:8}   '.format(datetime.date.today().strftime("%Y%m%d"))\
                    + '{:6}   '.format(datetime.datetime.utcnow().strftime("%H%M%S"))\
                    + '{:+.10e}'.format(chi2_ini/len(time_serie['dates']))\
                    + '\n'
                file_chains.write(format4chains)
                file_chains.close()

    try:
        del lnprobc
        del rstate
    except:
        pass

    # Create an archive
    path_code = cfgsetup.get('FullPaths', 'Code')
    path_event = cfgsetup.get('FullPaths', 'Event')
    path_arch = path_event + cfgsetup.get('RelativePaths', 'Archives')

    if os.path.exists(path_arch + cfgsetup.get('Controls', 'Archive')):
        shutil.rmtree(path_arch + cfgsetup.get('Controls', 'Archive'))

    dir = path_event + cfgsetup.get('RelativePaths', 'Chains')
    shutil.copytree(dir, path_arch + cfgsetup.get('Controls', 'Archive') + "/" + cfgsetup.get('RelativePaths', 'Chains'))

    dir = path_event + cfgsetup.get('RelativePaths', 'Data')
    shutil.copytree(dir, path_arch + cfgsetup.get('Controls', 'Archive') + "/" + cfgsetup.get('RelativePaths', 'Data'))

    file = path_event + "mulan.py"
    shutil.copyfile(file, path_arch + cfgsetup.get('Controls', 'Archive') + "/mulan.py")

    file = path_event + "observatories.ini"
    shutil.copyfile(file, path_arch + cfgsetup.get('Controls', 'Archive') + "/observatories.ini")

    file = path_event + "setup.ini"
    shutil.copyfile(file, path_arch + cfgsetup.get('Controls', 'Archive') + "/setup.ini")

    file = path_event + "advancedsetup.ini"
    shutil.copyfile(file, path_arch + cfgsetup.get('Controls', 'Archive') + "/advancedsetup.ini")

    dir = path_arch + cfgsetup.get('Controls', 'Archive') + "/" + cfgsetup.get('RelativePaths', 'Plots')
    os.makedirs(dir)

    try:
        shutil.rmtree(path_arch + cfgsetup.get('Controls', 'Archive') + "/" + cfgsetup.get('Controls', 'Archive') + ".zip")
    except:
        UserWarning("A ZIP file already exits. Archive not created.")

    filename = path_arch + cfgsetup.get('Controls', 'Archive')
    shutil.make_archive(filename, 'zip', filename)
    shutil.rmtree(filename)

    text = "Create archive {0:}".format(cfgsetup.get('Controls', 'Archive'))
    communicate(cfgsetup, 3, text, opts=False, prefix=False, newline=False, tab=True)

    # if os.path.exists(path_arch + cfgsetup.get('Controls', 'Archive') + '.zip'):
    #     os.remove(path_arch + cfgsetup.get('Controls', 'Archive') + '.zip')

    # shutil.move(filename + '.zip', path_arch)

    # Free memory
    try:
        del id_model, param_model, cond, id, accrate, key_list
        del value, chain_lengh, sampler, ndim, nwalkers
        del fitted_param, nuisance, result, params, grid_params
        del nb_params_grid
    except:
        pass

    if stop==1:
        sys.exit("\nProcess stopped by the user.\n")

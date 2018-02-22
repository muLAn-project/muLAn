# -*- coding: utf-8 -*-
"""Create a figure from muLAn outputs of publication quality"""

# Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan
# Distributed under the terms of the MIT license
#
# This module is part of software:
#       muLAn: gravitational MICROlensing Analysis code
#       https://github.com/muLAn-project/muLAn

import ConfigParser as cp
import glob
import numpy as np
from copy import copy
import muLAn.mulan as mulan
import muLAn.packages.general_tools as gtools
import matplotlib.pyplot as plt


class figure():
    """Create a figure from muLAn outputs
        
        Calling muLAnFig
        ================
        fig = muLAnFig.figure(figsize=(10,6), labelposx=0.8,
            labelposy=[0.6, 0.9], labelsize=10)
        
        Parameters
        ----------
        figsize : tuple, optional
            Figure size (x, y). Default is: figsize=(10,6).
        labelposx: float, optional
            Horizontal position of the labels of the observatories
            names. Possible values are 0 (left) to 1 (right).
            Default is: labelposx=0.8.
        labelposy: sequence of float, optional
            A 2-length sequence [ymin, ymax], which defines the
            vertical range where the labels of the observatories names
            will be displayed. Possible values are 0 (down) to 1 (up).
            Default is: labelposy=[0.6, 0.9].
        labelsize=10: float, optional
            Font size of the labels of the observatories names.
            Default is: labelsize=10
        
        User methods
        ------------
        addinset_caustics : see below
        addinset_lightcurve : see below
        save : see below
        show : see below

        Returns
        -------
        out : figure
            Using show will open a matplotlib.pylab interactive plot.
            Using save will store the figure in the disk. NB: do not
            use show before save, the saved figure will be empty.

        Examples
        --------
        Examples in automatic mode (called from the EVENT/ working directory):
        
        >>> fig = muLAnFig.figure(figsize=(10,6), labelposx=0.83,
                labelposy=[0.54, 0.94], labelsize=9)
        >>> fig.plot(trange=[7100, 7200], lcrange=[12, 16.8],
                resrange=[-0.38, 0.38])
        >>> fig.addinset_caustics([0.17, 0.64, 0.2, 0.3], xrange=[-1.75, -1.55],
                yrange=[-0.12, 0.13])
        >>> fig.addinset_lightcurve([0.15, 0.65, 0.3, 0.3], trange=[7144, 7148],
                lcrange=[12, 15.8])
        >>> fig.save('Plots/Figure.pdf')
        >>> fig.show()
        
        >>> fig = muLAnFig.figure()
        >>> fig.plot()
        >>> fig.addinset_caustics([0.2, 0.7, 0.2, 0.2])
        >>> fig.addinset_lightcurve([0.2, 0.4, 0.2, 0.2])
        >>> fig.show()
        
        Example in manual mode:
        
        >>> fig = muLAnFig.figure(labelposx=0.83, labelposy=[0.54, 0.94])
        >>> fig.plot(data=[('data1.dat', '#000000', 'Tel1'),
            ('data2.dat', '#FF00FF', 'Tel2')], lctraj=[('EARTH.dat', 'black')],
            trange=[7100, 7200], lcrange=[12, 16.8], resrange=[-0.38, 0.38])
        >>> fig.addinset_caustics([0.17, 0.64, 0.2, 0.3], caus=[('caustic.dat', 'red')])
        >>> fig.save('Plots/Figure.pdf')
        """
    def __init__(self, figsize=(10,6), labelposx=0.8, labelposy=[0.6, 0.9], labelsize=10):
        """Inititalize figure layout and search for muLAn output files"""
        self._labelposx = labelposx
        self._labelposy = labelposy
        self._labelsize = labelsize
        self._lccompt = 1
        self._causcompt = 1
        try:
            print " Searching for muLAn outputs...",
            self._getconfig()
            print "\033[1m\033[32mfound\033[0m"
        except:
            self._cfgsetup = None
            self._cfgobs = None
            print "\033[1m\033[35mnot found\033[0m (may lead to an error in non-manual mode)"
        try:
            print " Searching for muLAn best-fit parameters file...",
            self._getbestfitparams()
            print "\033[1m\033[32mfound\033[0m"
        except:
            print "\033[1m\033[35mnot found\033[0m (may lead to an error in non-manual mode)"
        # figure layout
        plt.close('all')
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        self.fig, (self._LC, self._RES) = plt.subplots(2, sharex=True, figsize=figsize, gridspec_kw={'height_ratios':[3, 1]})
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97, wspace=None, hspace=0.)
        
    def plot(self, data=None, lctraj=None, trange=None, lcrange=None, resrange=None):
        """Create main figure pannel
            
            Parameters
            ----------
            data: sequence of tuples, optional
                A n-length sequence [..., ('data_i.dat', 'color_i', 'label_i'), ...],
                where 'data_i.dat' is the ith data file, 'color_i' its color (name or
                hexadecimal code), and 'label_i' the name to be displayed on the
                figure labels. Default is: use muLAn outputs (data=None).
            lctraj: sequence of tuples, optional
                A n-length sequence [..., ('lctraj_i.dat', 'color_i'), ...], where
                'lctraj_i.dat' is the ith trajectory + light curve file, and
                'color_i' its color (name or hexadecimal code).
                Default is: use muLAn outputs (lctraj=None).
            trange: sequence of float, optional
                A 2-length sequence [tmin, tmax], which defines
                the range of dates to be plotted in the main plot.
                Default is: automatic range (trange=None)
            lcrange: sequence of float, optional
                A 2-length sequence [magmin, magmax], which defines
                the range of magnitudes to be plotted in the main plot.
                Default is: automatic range (lcrange=None)
            resrange: sequence of float, optional
                A 2-length sequence [magmin, magmax], which defines
                the range of magnitudes to be plotted in the residuals plot.
                Default is: automatic range (resrange=None)
            """
        print "\033[1m Creating main figure layout...\033[0m"
        # test whether to select outputs from muLAn's .ini files
        if (data == None) and (self._cfgobs != None):
            data = self._getdata()
        if (lctraj == None) and (self._cfgobs != None):
            lctraj = self._getlctraj()
        self.data = data
        self.lctraj = lctraj
        if trange: self._RES.set_xlim([trange[0], trange[1]])
        if resrange: self._RES.set_ylim([resrange[0], resrange[1]])
        if lcrange: self._LC.set_ylim([lcrange[0], lcrange[1]])
        self._LC.invert_yaxis()
        lwidth = 1.
        fontsize = 12
        plt.setp(self._LC.spines.values(), linewidth=lwidth)
        plt.setp(self._RES.spines.values(), linewidth=lwidth)
        self._LC.tick_params(labelsize=fontsize, width=0.8, direction='in', length=8)
        self._RES.tick_params(labelsize=fontsize, width=0.8, direction='in', length=8)
        self._RES.set_xlabel(r'HJD - 2,450,000', size=fontsize)
        self._LC.set_ylabel(r'Magnitude', size=fontsize)
        self._RES.set_ylabel(r'Residuals', size=fontsize)
        # plot theoretical light curves
        for lctraji, color in self.lctraj:
            print "   Reading theoretical light curve file:\033[3m", lctraji, "\033[0m"
            hjdr, amp, magr, xt, yt = np.loadtxt(lctraji, unpack=True)
            hjd, mag = self._optimizemc(hjdr, magr)
            self._LC.plot(hjd, mag, color=color, linewidth=1)
        # observatory names label positions
        if not self._labelposx: x = 0.7
        else: x = self._labelposx
        if not self._labelposy: ymin, ymax = 0.5, 0.7
        else: ymin, ymax = self._labelposy
        y = iter(np.linspace(ymin, ymax, len(self.data), endpoint=False))
        # read and plot data and residuals
        for datai, color, obsname in self.data:
            print "   Reading data file:\033[3m", datai, "\033[0m"
            ID, hjd, mag, errmag, resmag, amp, erramp, resamp, bkg, seeing, xs, ys, chi, toto, tata = np.loadtxt(datai, unpack=True)
            self._LC.errorbar(hjd, mag, errmag, fmt='o', color=color, markersize=4, alpha=0.4, linewidth=1)
            self._RES.plot((np.min(hjd), np.max(hjd)), (0., 0.), 'k-', linewidth=0.4)
            self._RES.errorbar(hjd, resmag, errmag, fmt='o', color=color, markersize=4, alpha=0.4, linewidth=1)
            # display observatory names
            self._LC.annotate(r'$\bullet$ ' + obsname, xy=(x, y.next()), xycoords='figure fraction', color=color, fontsize=self._labelsize)

    def addinset_lightcurve(self, layout, trange=None, lcrange=None):
        """Add a light curve zoom pannel to the plot
            
            Parameters
            ----------
            layout: sequence of float
                A 4-length sequence [left, bottom, width, height], which
                defines the position and shape of the inset.
            trange: sequence of float, optional
                A 2-length sequence [tmin, tmax], which defines
                the range of dates to be plotted in the inset plot.
                Default is: automatic range (trange=None)
            lcrange: sequence of float, optional
                A 2-length sequence [magmin, magmax], which defines
                the range of magnitudes to be plotted in the inset plot.
                Default is: automatic range (lcrange=None)
            """
        print "\033[1m Creating light curve inset " + str(self._lccompt) + "...\033[0m"
        self._lccompt += 1
        # pannel layout
        ZLC = self.fig.add_axes(layout)
        if trange: ZLC.set_xlim(trange)
        if lcrange: ZLC.set_ylim(lcrange)
        ZLC.invert_yaxis()
        # read and plot observed data and residuals
        for datai, color, obsname in self.data:
            print "   Reading data file:\033[3m", datai, "\033[0m"
            ID, hjd, mag, errmag, resmag, amp, erramp, resamp, bkg, seeing, xs, ys, chi, toto, tata = np.loadtxt(datai, unpack=True)
            ZLC.errorbar(hjd, mag, errmag, fmt='o', color=color, markersize=4, alpha=0.4, linewidth=1)
        # plot theoretical light curves
        for lctraji, color in self.lctraj:
            print "   Reading theoretical light curve file:\033[3m", lctraji, "\033[0m"
            hjdr, amp, magr, xt, yt = np.loadtxt(lctraji, unpack=True)
            hjd, mag = self._optimizemc(hjdr, magr)
#            hjd, amp, mag, xt, yt = np.loadtxt(lctraji, unpack=True)
            ZLC.plot(hjd, mag, color=color, linewidth=1)

    def addinset_caustics(self, layout, caus=None, xrange=None, yrange=None):
        """Add a caustic pannel to the plot
            
            Parameters
            ----------
            layout: sequence of float
                A 4-length sequence [left, bottom, width, height], which
                defines the position and shape of the inset.
            caus: sequence of tuples, optional
                A n-length sequence [..., ('caus_i.dat', 'color_i'), ...], where
                'caus_i.dat' is the ith caustic file, and 'color_i' its color (name or
                hexadecimal code). Default is: use muLAn outputs (caus=None).
            xrange: sequence of float, optional
                A 2-length sequence [xmin, xmax], which defines
                the horizontal range for the caustic plot.
                Default is: automatic range (xange=None)
            yrange: sequence of float, optional
                A 2-length sequence [ymin, ymax], which defines
                the vertical range for the caustic plot.
                Default is: automatic range (yrange=None)
            """
        print "\033[1m Creating caustics inset "+ str(self._causcompt) + "...\033[0m"
        self._causcompt += 1
        # pannel layout
        CAU = self.fig.add_axes(layout)
        if not (xrange and yrange): CAU.set_aspect('equal')
        if xrange: CAU.set_xlim(xrange)
        if yrange: CAU.set_ylim(yrange)
        # plot trajectories
        for lctraji, color in self.lctraj:
            print "   Reading trajectory file:\033[3m", lctraji, "\033[0m"
            hjd, amp, mag, xt, yt = np.loadtxt(lctraji, unpack=True)
            CAU.plot(xt, yt, color=color, linewidth=1)



        # Load caustics
        fname = "{:s}CAUSTIC.dat".format(self._pathoutputs)
        fcaustics = np.loadtxt(fname, unpack=False, dtype=np.float64)
        n_caus = fcaustics.shape[1] / 2

        # Load caustic times
        if n_caus > 1:
            fname = "{:s}CAUSTIC.dat".format(self._pathoutputs)
            fcaustics = open(fname, 'r')
            for line in fcaustics: break
            fcaustics.close()
            times = line.replace(" ", "").replace("x", "").replace("y", "").replace("#", "").replace("(", "").replace("\n", "").split(")")
            times = np.unique(times)
            try:
                times[0] = 't0 = {:.6f}'.format(self._bf['t0'])
            except AttributeError:
                times[0] = "t0"
        else:
            try:
                times = np.atleast_1d(['t0 = {:.6f}'.format(self._bf['t0'])])
            except AttributeError:
                times = np.atleast_1d(['t0'])

        # Plot caustics
        if caus == None:
            color_caus = ['red', 'Orange', 'SeaGreen', 'LightSeaGreen', 'CornflowerBlue', 'DarkViolet']
        else:
            color_caus = [color for cau, color in caus]

        for i in range(n_caus):
            # print "   Plotting caustic " + str(i + 1) + "..."
            print "   Plotting caustic:\033[3m", times[i], "\033[0m"
            xc = fcaustics.T[2*i]
            yc = fcaustics.T[2*i + 1]
            CAU.scatter(xc, yc, marker='.', c=color_caus[i], s=0.1)
            color_caus = np.roll(color_caus, -1)

    def save(self, figname):
        """Save figure"""
        print "\033[1m Saving figure: \033[3m" + figname + "...\033[0m"
        plt.savefig(figname)

    def show(self):
        """Show figure"""
        plt.show()
    
    def _optimizemc(self, x, y, err=0.001):
        """Optimize the sampling of the input curve"""
        print "   ... optimizing sampling:",
        N = len(x)
        ts = np.zeros(N, dtype=np.float_)
        As = np.zeros(N, dtype=np.float_)
        # Algorithm
        n = 0
        i = 0
        while n < N:
            cond = True
            As[n] = y[i]
            ts[n] = x[i]
            n += 1
            p = 2
            while p <= N-1-i: # 2≤p
                if np.logical_not(cond):
                    break
                for k in np.arange(p-1)+1: # 1≤k≤p-1
                    Alin = (x[i+k] - x[i]) * (y[i+p] - y[i]) / (x[i+p] - x[i]) + y[i]
                    cond = np.abs(Alin - y[i+k]) <= err
                    if np.logical_not(cond):
                        i = i+p-1
                        break
                p += 1
            if (p == N-i): break
        ts[n-1] = x[i]
        As[n-1] = y[i]
        ts[n] = x[N-1]
        As[n] = y[N-1]
        xopt = copy(ts[0:n+1])
        yopt = copy(As[0:n+1])
        # verbose
        print "\033[3mkeeping " + str(n + 1) + " points out of " + str(N) + "\033[0m"
        return xopt, yopt

    def _getdata(self):
        """Get usefull data file names from muLAn's .ini files"""

        # Extract required information
        data = list()
        for i in xrange(len(self._observatories)):
            table = [a.strip() for a in self._cfgobs.get('ObservatoriesDetails', self._observatories[i]).split(',')]
            fname = "{:s}{:s}.dat".format(self._pathoutputs, self._observatories[i].upper())
            data.append((fname, "#" + table[1], table[0]))

        return data

    def _getlctraj(self):
        """Get usefull light curve/trajetory file names from muLAn's .ini files"""

        # Extract required information
        traj = list()
        colors = np.array(['black', 'blue', 'orange'])
        for i in xrange(len(self._observatories)):
            table = [a.strip() for a in self._cfgobs.get('ObservatoriesDetails', self._observatories[i]).split(',')]
            traj.append(table[4])
        traj = np.unique(traj)

        traj_final = list()
        for a in traj:
            fname = "{:s}{:s}.dat".format(self._pathoutputs, a.upper())
            traj_final.append((fname, colors[0]))
            colors = np.roll(colors, -1)

        return traj_final

    def _getconfig(self):
        """Load the configuration files *.ini."""

        # Path of the event
        path_event = mulan.getpath_event()

        # Configuration files
        fname_setup = "{:s}setup.ini".format(path_event)
        fname_advanced = "{:s}advancedsetup.ini".format(path_event)
        fname_obs = "{:s}observatories.ini".format(path_event)

        # Load configuration files
        cfgsetup = cp.SafeConfigParser()
        cfgsetup.read([fname_setup, fname_advanced])
        cfgobs = cp.SafeConfigParser()
        cfgobs.read(path_event + 'observatories.ini')

#        text = "Load parameter files..."
#        gtools.communicate(cfgsetup, 1, text, opts=[gtools.printoption.level0], prefix=True, newline=True)

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

        self._cfgsetup = cfgsetup
        self._cfgobs = cfgobs

        # Data file names
        data2find = np.array(self._cfgobs.options('ObservatoriesDetails'))

        # Cross-match with data to use
        data2use = np.array(self._cfgsetup.options('Observatories'))
        data2use = np.array([data2find[i] for i in xrange(len(data2find)) if np.any(data2use == data2find[i])])

        # Cross-match with existing files
        path = self._cfgsetup.get('FullPaths', 'Event') + self._cfgsetup.get('RelativePaths', 'Outputs')
        data_filenames = glob.glob(path + '*')
        observatories = [a.split('/')[-1] for a in data_filenames]
        data_filenames = [data_filenames[i] for i in xrange(len(data_filenames)) if
                          np.any(data2use == observatories[i].rpartition('.')[0].lower())]
        observatories = [ob.rpartition('.')[0].lower() for ob in observatories if
                         np.any(data2use == ob.rpartition('.')[0].lower())]
        self._observatories = observatories
        self._pathoutputs = path

    def _getbestfitparams(self):
        """Load the values of the best fit.""" 
    
        fname = "{:s}Results.txt".format(self._pathoutputs)
        file_res = open(fname, 'r')
        res = ""
        for line in file_res:
            res = res + line
        file_res.close()

        res = res.split("Best-fitting parameters")[1]
        res = res.split("Site")[0]
        res = res.replace("\n", "").split(" ")
        res = [a for a in res if a != ""]
        res = [a for a in res if a != "="]
        res = np.reshape(res, (len(res)/2, 2))

        bf = dict()
        [bf.update({res.T[0][i]: float(res.T[1][i])}) for i in range(res.shape[0])]
        self._bf = bf

if __name__ == '__main__':
    help(muLAnFig)


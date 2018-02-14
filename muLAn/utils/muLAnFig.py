# -*- coding: utf-8 -*-
"""muLAnFig: a tool to create a figure from muLAn outputs"""

# Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan
# Distributed under the terms of the MIT license
#
# This module is part of software:
#       muLAn: gravitational MICROlensing Analysis code
#       https://github.com/muLAn-project/muLAn

import numpy as np
import matplotlib.pyplot as plt


class figure():
    """Create a figure from muLAn outputs
        
        Calling muLAnFig
        ================
        mlf = muLAnFig(data=None, lctraj=None, caus=None, trange=None, lcrange=None, resrange=None, labelpos=None, labelsize=10, figsize=(10,6))
        
        Parameters
        ----------
        (coming soon)

        Examples
        --------
        (coming soon)
        """
    def __init__(self, figsize=(10,6), labelposx=None, labelposy=None, labelsize=10):
        """Create figure layout"""
        self._labelposx = labelposx
        self._labelposy = labelposy
        self._labelsize = labelsize
        # figure layout
        plt.close('all')
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        self.fig, (self._LC, self._RES) = plt.subplots(2, sharex=True, figsize=figsize, gridspec_kw={'height_ratios':[3, 1]})
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.97, top=0.97, wspace=None, hspace=0.)
        
    def plot(self, data=None, lctraj=None, caus=None, trange=None, lcrange=None, resrange=None):
        """Create main figure pannel"""
        print "\033[1m Creating main figure layout...\033[0m"
        # test whether to select outputs from muLAn's .ini files
        if not data:
            data = _getdata()
        if not lctraj:
            traj = _getlctraj()
        if not caus:
            caus = _getcaus()
        self.data = data
        self.lctraj = lctraj
        self.caus = caus
        self._lccompt = 1
        self._causcompt = 1
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
            hjd, amp, mag, xt, yt = np.loadtxt(lctraji, unpack=True)
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
            self._LC.errorbar(hjd, mag, errmag, fmt='o', color=color, markersize=4, alpha=0.3, linewidth=1)
            self._RES.plot((np.min(hjd), np.max(hjd)), (0., 0.), 'k-', linewidth=0.4)
            self._RES.errorbar(hjd, resmag, errmag, fmt='o', color=color, markersize=4, alpha=0.3, linewidth=1)
            # display observatory names
            self._LC.annotate(r'$\bullet$ ' + obsname, xy=(x, y.next()), xycoords='figure fraction', color=color, fontsize=self._labelsize)

    def addinset_lightcurve(self, layout, trange=None, lcrange=None):
        """Add a light curve zoom pannel to the plot
            
            Parameters
            ----------
            (coming soon)
            
            Examples
            --------
            (coming soon)
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
            ZLC.errorbar(hjd, mag, errmag, fmt='o', color=color, markersize=4, alpha=0.3, linewidth=1)
        # plot theoretical light curves
        for lctraji, color in self.lctraj:
            print "   Reading theoretical light curve file:\033[3m", lctraji, "\033[0m"
            hjd, amp, mag, xt, yt = np.loadtxt(lctraji, unpack=True)
            ZLC.plot(hjd, mag, color=color, linewidth=1)

    def addinset_caustics(self, layout, caus=None, xrange=None, yrange=None):
        """Add a caustic pannel to the plot
            
            Parameters
            ----------
            (coming soon)

            Examples
            --------
            (coming soon)
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
        # plot caustics
        for causi, color in self.caus:
            print "   Reading caustic file:\033[3m", causi, "\033[0m"
            xc, yc = np.loadtxt(causi, unpack=True)
            CAU.scatter(xc, yc, marker='.', c='red', s=0.1)

    def save(self, figname):
        """Save figure"""
        print "\033[1m Saving figure: \033[3m" + figname + "...\033[0m"
        plt.savefig(figname)

    def show(self):
        """Show figure"""
        plt.show()

    def _getdata(self):
        """Get usefull data file names from muLAn's .ini files"""
        # merci d'utiliser le format: list( tuple('file path', 'color', 'telescope'), ... )
        # e.g.: [('Outputs/OGLE-I.dat', '#FFFF33', 'OGLE'), ...]
        data = list()
        return data

    def _getlctraj(self):
        """Get usefull light curve/trajetory file names from muLAn's .ini files"""
        # merci d'utiliser le format: list( tuple('file path', 'color'), ... )
        # e.g.: [('Outputs/EARTH.dat', '#FFFF33'), ...]
        traj = list()
        return traj

    def _getcaus(self):
        """Get usefull caustics file names from muLAn's .ini files"""
        # merci d'utiliser le format: list( tuple('file path', 'color'), ... )
        # e.g.: [('Outputs/caustics_t1.dat', '#FFFF33'), ...]
        caus = list()
        return caus


if __name__ == '__main__':
    help(muLAnFig)


# -*- coding: utf-8 -*-
"""muLAnCleanData: an interactive tool to select data points from muLAn input files for cleaning purpose"""

# Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan
# Distributed under the terms of the MIT license
#
# This module is part of software:
#       muLAn: gravitational MICROlensing Analysis code
#       https://github.com/muLAn-project/muLAn
#
# This module:
#       Copyright (c) 2018 Yiannis Tsapras
#       Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan

import numpy as np
import matplotlib.pyplot as plt
from sys import argv, exit


class PointBrowser():
    """Browse point on screen"""
    def __init__(self, fig, ax, line, idx, hjd, mag):
        self.fig = fig
        self.ax = ax
        self.line = line
        self.idx = idx
        self.hjd = hjd
        self.mag = mag
        self.lastind = 0
        self.text = self.ax.text(0.05,0.95,'selected: none',
                            transform=self.ax.transAxes, va='top')
        self.selected, = self.ax.plot([self.hjd[0]], [self.mag[0]], 'k.', ms=12, alpha=0.5,\
                                 color='yellow', visible=False)
        self.output = []
    
    def onpress(self, event):
        if self.lastind is None:
            return
        if event.key not in ('n','p','a','r','t','left','right'):
            return
        if (event.key == 'n' or event.key == 'right'):
            inc = 1
        if (event.key == 'p' or event.key == 'left'):
            inc = -1
        if (event.key == 'a'):
            if (self.idx[self.lastind]) not in self.output:
                self.output.append(self.idx[self.lastind])
            return
        if (event.key == 'r'):
            if (self.idx[self.lastind]) in self.output:
                self.output.remove(self.idx[self.lastind])
            return
        if event.key == 't':
            the_list = self.output
            the_list.sort()
            print (the_list)
            return
        self.lastind += inc
        self.lastind = np.clip(self.lastind, 0, len(hjd) - 1)
        self.update()
      
    def onpick(self, event):
        if event.artist != self.line:
            return True
        N = len(event.ind)
        if not N:
            return True
        # click position
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        # select point closest to the position clicked
        distances = np.hypot(x - self.hjd[event.ind], y - self.mag[event.ind])
        indmin = distances.argmin()
        dataind = event.ind[indmin]
        self.lastind = dataind
        self.update()
    
    def update(self):
        if self.lastind is None:
            return
        dataind = self.lastind
        self.ax.plot(self.hjd, self.mag, 'k.')
        self.ax.plot(self.hjd[self.output], self.mag[self.output], 'r.')
        self.selected.set_visible(True)
        self.selected.set_data(self.hjd[dataind], self.mag[dataind])
        self.text.set_text("selected: %s" % str(self.idx[dataind]))
        self.fig.canvas.draw()

def cleandata(datafile):
    """muLAn data cleaning tool
        
        Click on a point to highlight it. Navigate through the points with the
        points using the 'n'/'right' and 'p'/'left' keys for next and previous
        points. Add a new point to the list by pressing 'a' or remove it by
        pressing 'r'. Print the list at any point by pressing 't'.
        Press 'q' to quit.
        """
    print "\033[3m "+ cleandata.__doc__ + "\033[0m"
    # read data file
    try:
        idx, hjd, mag, merr, fwhm, bkg = np.loadtxt(datafile, unpack=True)
        idx = idx.astype(int)
    except: raise IOError("unable to read %s" % datafile)
    # create interactive plot
    fig, ax = plt.subplots()
    line, = ax.plot(hjd, mag, 'k.', picker=5) # 5 points tolerance
    plt.ylim(max(mag) + 0.1, min(mag) - 0.1)
    browser = PointBrowser(fig, ax, line, idx, hjd, mag)
    fig.canvas.mpl_connect('pick_event', browser.onpick)
    fig.canvas.mpl_connect('key_press_event', browser.onpress)
    plt.show()

if __name__ == '__main__':
    help(muLAnCleanData)
    

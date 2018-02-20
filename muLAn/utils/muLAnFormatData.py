# -*-coding:Utf-8 -*
"""formatdata: a tool to convert data files in muLAn format"""

# Copyright (c) 2014-2018 Clément Ranc & Arnaud Cassan
# Distributed under the terms of the MIT license
#
# This module is part of software:
#       muLAn: gravitational MICROlensing Analysis code
#       https://github.com/muLAn-project/muLAn

import numpy as np
import os


def formatdata(infilename, outfilename, cols):
    """Reformat data files to be used by muLAn
        
        Calling formatdata
        ==================
        formatdata(infilename, outfilename, cols)
        
        Usage
        -----
        Enter in cols the list of columns description (i.e. keywords,
        see below) in the order they appear in the input file.
        
        Parameters
        ----------
        infilename: string
            Name of input data file.
        outfilename: string
            Name of output data file in muLAn format.
        cols: sequence of strings
            Mandatory keywords are:
                'hjd': Julian date or modified Julian date.
                'mag': magnitude.
                'errmag': error in magnitude.
            Optional keywords are:
                'seeing': seeing.
                'backg': background.
            For useless columns, use e.g. 'other'

        Examples
        --------
        >>> formatdata('data.dat', 'data_muLAn.dat',
                ['hjd', 'mag', 'errmag', 'seeing', 'backg'])
        
        >>> formatdata('data.dat', 'data_muLAn.dat',
                ['other', 'hjd', 'mag', 'errmag'])
        """
    # check mandatory keywords
    mandatkeys = ['hjd', 'mag', 'errmag']
    for key in mandatkeys:
        # check whether all mandatory keywords are present
        if key not in cols:
            raise ValueError("mandatory column missing: " + key)
        # check whether keywords appear only once
        if cols.count(key) > 1:
            raise ValueError("column appears more than once: " + key)
    # check if input file exists, and process it
    if not os.path.isfile(infilename):
        raise IOError("file '" + infilename + "' does not exist")
    # limit number of columns to read
    usecols = range(len(cols))
    # reading input data file
    print " Reading input data file:\033[3m", infilename, "\033[0m"
    dtype = {'names': tuple(cols), 'formats': tuple(['S50' for c in cols])}
    data = np.loadtxt(infilename, dtype=dtype, usecols=usecols, unpack=False)
    # re-order columns
    newfile = ''
    for i in range(len(data['hjd'])):
        # check whether date is in HJD or MHJD, and correct it
        mhjd = float(data['hjd'][i]) - 2450000.
        if mhjd > 0.:
            data['hjd'][i] = str(mhjd)
        # mandatory keywords
        newfile = newfile + repr(i + 1) + ' ' + data['hjd'][i] + ' ' + data['mag'][i] + ' ' + data['errmag'][i]
        # optional keywords
        if 'seeing' in cols:
            newfile = newfile + ' ' + data['seeing'][i]
        else:
            newfile = newfile + ' 0'
        if 'backg' in cols:
            newfile = newfile + ' ' + data['backg'][i] + '\n'
        else:
            newfile = newfile + ' 0\n'
    # create output data file in muLAn format
    print " Creating output data file in muAn format:\033[3m", outfilename, "\033[0m"
    outfile = open(outfilename, 'w')
    outfile.write(newfile)
    outfile.close()

if __name__ == '__main__':
    help(cleandata)


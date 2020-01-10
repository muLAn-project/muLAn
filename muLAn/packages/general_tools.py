# -*-coding:Utf-8 -*
import os
import sys
import numpy as np
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
# --------------------------------------------------------------------
# Fonctions
# --------------------------------------------------------------------
def communicate(cfg, verbose, text, opts=False, prefix=False, newline=False, tab=False):
    if cfg.getint('Modelling', 'Verbose') >= verbose:
        if prefix:
            text = "[muLAn] " + text
        if opts!=False:
            text2=''
            for a in opts:
                text2 = text2 + a
                text = text2 + text + '\033[0m'
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print(text)
        else:
            if tab:
                text = "    " + text
            if newline:
                text = "\n" + text
            print(text)
# --------------------------------------------------------------------
def add_topythonpath(path):
    if len(path)>0:
        sys.path.insert(0, path)
# --------------------------------------------------------------------
def check_slash(path):
    if len(path) > 0:
        if path[-1] != '/':
            path = '{:s}/'.format(path)
    return path
# --------------------------------------------------------------------
def unpack_options(cfgsetup, level0, level1):
    options = [a.strip() for a in cfgsetup.get(level0, level1).split(',')]
    del a, cfgsetup, level0, level1
    return options
# ----------------------------------------------------------------------
def str2bool(string):
    if (string == 'True') | (string == 'true') | (string == '1'):
        return True
    else:
        return False

# ----------------------------------------------------------------------
def print_params(list):
    text = ""
    if list[0]=='fit':
        a = np.abs(float(list[1]))
        b = np.abs(float(list[2]))
        c = float(list[3])
        text = "{:.6f} chosen within {:.6f} --> {:.6f}".format(c, c-a, c+b)
    if list[0]=='fix':
        c = float(list[3])
        text = "{:.6f} (fixed)".format(c)
    if list[0]=='gri':
        a = np.abs(float(list[1]))
        b = np.abs(float(list[2]))
        c = int(list[3])
        text = "on a grid from {:.6f} to {:.6f} ({:d} points)".format(a, b, c)
    return text

# ----------------------------------------------------------------------
def prompt(path, name):
    filename = path + name + '.zip'
    if os.path.exists(filename):
        text = 'The archive:\n   ' + filename + '\ndoes exist. \033[44m\033[97m\033[1mDo you want to overwrite it? [y/N]\033[0m '
        x = raw_input(text)
        if (x == 'N') | (x == 'n') | (x == ''):
            text = '\033[44m\033[97m\033[1mPlease enter a new filename (without extension): [ARCHIVE]\033[0m '
            y = raw_input(text)
            if (y == ''): y = 'ARCHIVE'
            name = prompt(path, y)
            return name
        else:
            return name
    else:
        return name

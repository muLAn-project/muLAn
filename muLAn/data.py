# -*-coding:Utf-8 -*
# ====================================================================
# Packages
# ====================================================================
#import configparser as cp
#import copy
#import glob
#import muLAn
#import muLAn.packages.general_tools as gtools
import numpy as np
#import os
import pandas as pd
#import sys

class Data:
    """Class to store and modify microlensing data.

    This class is not fully developed yet. Its purpose is to store all the
    relevant information about each observation that is used to fit a light
    curve.

    Args:
        data (dict or :obj:`pandas.DataFrame`): observations
    
    Attributes:

    """

    def __init__(self, data, **kwargs):

        if isinstance(data, dict):
            self.table = pd.DataFrame.from_dict(data).astype({
                'amp': np.dtype('f8'), 
                'fs': np.dtype('f8'),
                'fb': np.dtype('f8')})
        elif isinstance(data, pd.DataFrame):
            self.table = data


if (__name__ == "__main__"):
    pass


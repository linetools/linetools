""" Utilities for isgm
 Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb
import numpy as np
import warnings

from linetools import utils as ltu
from linetools.isgm.abssystem import GenericAbsSystem

def abssys_from_json(filename):
    """
    Parameters
    ----------
    filename

    Returns
    -------
    abs_sys : AbsSystem

    """
    # Load JSON file to determine type
    adict = ltu.loadjson(filename)
    if 'class' in adict.keys():
        if adict['class'] == 'MgIISystem':
            from pyigm.abssys.igmsys import MgIISystem
            abs_sys = MgIISystem.from_dict(adict)
        else:
            warnings.warn("Unknown or uncoded class: {:s}.\nMaking a Generic one".format(adict['class']))
            abs_sys = GenericAbsSystem.from_dict(adict)
    else:
        abs_sys = GenericAbsSystem.from_dict(adict)

    # Return
    return abs_sys

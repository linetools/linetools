"""
Module for utilites related to analysis of lines
  -- Intended to be methods generic to emission and absorption (e.g. EW)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy import constants as const


#from xastropy.xutils import xdebug as xdb

# EW 
def box_ew(spec):
    """  Boxcar EW calculation
    Observer frame, not rest-frame

    Parameters
    ----------
    spec -- Tuple of (wave, fx, sig)

    Returns:
    ----------
      EW, sigEW : EW and error in observer frame
    Note: Tested in test_absline_anly
    """
    # Grab
    wv,fx,sig = spec

    # Cut spectrum
    # dwv
    dwv = wv - np.roll(wv,1)
    dwv[0] = dwv[1]


    # Simple boxcar
    EW = np.sum( dwv * (1. - fx) ) 
    varEW = np.sum( dwv**2 * sig**2 )
    sigEW = np.sqrt(varEW) 

    # Return
    return EW, sigEW

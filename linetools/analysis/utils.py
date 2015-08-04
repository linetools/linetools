"""
Module for utilites related to analysis of lines
  -- Intended to be methods generic to emission and absorption (e.g. EW)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy import constants as const
from astropy.modeling import models, fitting

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

#Gaussian EW
def gaussian_ew(spec):
    """  EW calculation using Gaussian fit
    Observer frame, not rest-frame
    wvlim must be set!
    spec must be set!

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

    # dwv
    dwv = wv - np.roll(wv,1)
    dwv[0] = dwv[1]
    
    # Fit the data using a Gaussian
    # Initial guesses
    amp_init = np.mean(fx)/2. #half the mean flux
    stddev_init = 3*np.mean(dwv) #3 pixels
    mean_init = np.mean(wv) #half wave range
    # Model initialization
    g_init = models.GaussianAbsorption1D(amplitude=amp_init.value, mean=mean_init.value, stddev=stddev_init.value)
    # Fitting initialization
    fit_g = fitting.LevMarLSQFitter()
    # Actual fit
    g = fit_g(g_init, wv, fx, weights=1./sig)

    # Area under curve of Gaussian is [amplitude*stdev*sqrt(2*pi)]
    EW = g.amplitude.value * g.stddev.value * np.sqrt(2 * np.pi) #unitless
    EW = EW * wv.unit #add the same unit as wv
    #error missing; maybe sum up the residuals in cuadrature?
    sigEW = np.nan

    #Return
    return EW, sigEW

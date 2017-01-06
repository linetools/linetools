""" Line analysis tools

These are intended to be methods generic to emission and absorption
(e.g. Equivalent width)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy.modeling import models, fitting


def box_ew(spec):
    """  Boxcar EW calculation

    Observer frame, not rest-frame

    Parameters
    ----------
    spec : Tuple of (wave, fx, sig)

    Returns
    -------
    EW, sigEW : EW and error in observer frame
    """
    # Note: Tested in test_absline_anly

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


def gaussian_ew(spec, ltype, initial_guesses=None):
    """  EW calculation using Gaussian fit

    Observer frame, not rest-frame. wvlim and spec must be set!

    Parameters
    ----------
    spec : Tuple of (wave, fx, sig)
    ltype : string
      whether this is for absorption or emission line (see SpectralLine Class)
    initial_guesses, optional : Tuple of (amplitude, mean, stddev) 
      Initial guesses of the Gaussian fit (unitless)

    Returns
    -------
    EW, sigEW : EW and error in observer frame
    
    """
    # Note: Tested in test_absline_anly

    # Grab
    wv,fx,sig = spec

    # dwv
    dwv = wv - np.roll(wv,1)
    dwv[0] = dwv[1]

    # Initial guesses of the Gaussian fit
    if initial_guesses is None:
        amp_init = np.mean(fx).value/2. #half the mean flux
        stddev_init = 3*np.mean(dwv).value #3 pixels
        mean_init = np.mean(wv).value #half wave range
    elif len(initial_guesses)==3:
        amp_init = initial_guesses[0]
        mean_init = initial_guesses[1] 
        stddev_init = initial_guesses[2]
        #check whether these values are sensible
        if (mean_init < np.min(wv.value)) or (mean_init > np.max(wv.value)):
             raise ValueError('gaussian_ew: The initial guess for Gaussian mean is not sensible; check it!')
        if (amp_init < 0):
             raise ValueError('gaussian_ew: The initial guess for Gaussian amplitude is not sensible; check it!')
        if (stddev_init < 0):
             raise ValueError('gaussian_ew: The initial guess for Gaussian stddev is not sensible; check it!')
    else:
        raise ValueError('gaussian_ew: Format of the initial_guesses is incorrect')

    # Model initialization
    if ltype == 'Abs':
        g_init = models.GaussianAbsorption1D(amplitude=amp_init, mean=mean_init, stddev=stddev_init) # This model does not support units
    elif ltype == 'Emiss':
        g_init = models.Gaussian1D(amplitude=amp_init, mean=mean_init, stddev=stddev_init) # This model does not support units
    else:
        raise ValueError("gaussian_ew: ltype has to be either 'Abs' or 'Emiss'")    
    
    # Fitting algorithm initialization
    fit_g = fitting.LevMarLSQFitter()
    # Use only good values (i.e. with meaningful errors)
    cond = (sig > 0.) & (np.isfinite(sig))
    # Actual fit
    g = fit_g(g_init, wv[cond], fx[cond], weights=1./sig[cond])

    #Check whether the final fit is sensible
    fit_info = fit_g.fit_info
    if fit_info['param_cov'] is None:
        raise ValueError('gaussian_ew: The fit is not sensible! Check initial_guesses')

    # Area under curve of Gaussian is [amplitude*stddev*sqrt(2*pi)]
    EW = g.amplitude.value * g.stddev.value * np.sqrt(2 * np.pi) #unitless
    EW = EW * wv.unit #add the same unit as wv
    
    #error estimation
    cov = fit_g.fit_info['param_cov'] #covariance matrix
    x = g.parameters[0] # amplitude
    y = g.parameters[2] # stddev
    sigEW = EW * np.sqrt(cov[0,0] / x**2 + cov[2,2] / y**2 + 2 * cov[0,2] / (x*y))

    return EW, sigEW

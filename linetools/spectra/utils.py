""" Methods related to spectra
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import pdb

from astropy import units as u

from .xspectrum1d import XSpectrum1D

def splice_two(spec1, spec2, wvmx=None, scale=1., chk_units=True):
    """ Combine two overlapping spectra.

    Parameters
    ----------
    spec1 : XSpectrum1D
    spec2 : XSpectrum1D
      The overlapping spectrum. It should cover wavelengths
      *longer* than the original spectrum.
    wvmx : Quantity, optional
      Wavelength to begin splicing *after*
    scale : float, optional
      Scale factor for flux and error array.
      Mainly for convenience of plotting
    chk_units : bool, optional
      Check units of input spectra

    Returns
    -------
    spec3 : XSpectrum1D
      A copy of the spliced spectrum.
    """
    # Checks
    if chk_units:
        for key,item in spec1.units.items():
            assert item == spec2.units[key]
    # Begin splicing after the end of the internal spectrum
    if wvmx is None:
        wvmx = spec1.wvmax
    #
    gdp = np.where(spec2.wavelength > wvmx)[0]
    # Concatenate
    new_wv = np.concatenate((spec1.wavelength.value,
                             spec2.wavelength.value[gdp]))
    uwave = u.Quantity(new_wv, unit=spec1.units['wave'])
    new_fx = np.concatenate((spec1.flux.value,
                             spec2.flux.value[gdp] * scale))
    # Error
    if spec1.sig_is_set:
        new_sig = np.concatenate((spec1.sig, spec2.sig[gdp] * scale))
    else:
        new_sig = None

    # Continuum
    if spec1.co_is_set:
        new_co = np.concatenate((spec1.co, spec2.co[gdp] * scale))
    else:
        new_co = None

    # Generate
    spec3 = XSpectrum1D.from_tuple(
        (uwave, u.Quantity(new_fx), new_sig, new_co), meta=spec1.meta.copy())
    # Return
    return spec3

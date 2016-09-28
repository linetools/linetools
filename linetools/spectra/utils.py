""" Methods related to spectra
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import json
import pdb

from astropy import units as u
from astropy import constants as const

from linetools import utils as liu


def meta_to_disk(in_meta):
    """ Polish up the meta dict for I/O

    Parameters
    ----------
    in_meta : dict

    Returns
    -------
    meta_out : str

    """
    meta = in_meta.copy()
    # Headers
    for kk,header in enumerate(in_meta['headers']):
        if header is None:
            meta['headers'][kk] = str('none')
        else:
            try:
                meta['headers'][kk] = header.tostring()
            except AttributeError:
                if not isinstance(header, basestring):
                    raise ValueError("Bad format in header")
                meta['headers'][kk] = header
    # Clean up the dict
    d = liu.jsonify(meta)
    return json.dumps(d)


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
    from .xspectrum1d import XSpectrum1D

    # Checks
    if chk_units:
        for key,item in spec1.units.items():
            assert item == spec2.units[key]
    if spec2.wvmax < spec1.wvmax:
        raise IOError("The second input spectrum does not cover longer wavelengths.")
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


def rebin_to_rest(spec, zarr, dv, debug=True):
    """ Shuffle an XSpectrum1D dataset to an array of
    observed wavelengths and rebin to dv pixels.
    Note: This works on the unmasked data array and returns
    unmasked spectra

    Parameters
    ----------
    spec : XSpectrum1D
    zarr : nd.array
      Array of redshifts
    dv : Quantity
      Velocity width of the new pixels

    Returns
    -------
    new_spec : XSpectrum1D
      Not masked

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D

    # Error checking
    if spec.nspec <= 1:
        raise IOError("Use spec.rebin instead")
    if spec.nspec != len(zarr):
        raise IOError("Input redshift array must have same dimension as nspec")
    # Generate final wave array
    dlnlamb = np.log(1+dv/const.c)
    z2d = np.outer(zarr, np.ones(spec.totpix))
    wvmin, wvmax = (np.min(spec.data['wave']/(1+z2d))*spec.units['wave'],
                    np.max(spec.data['wave']/(1+z2d))*spec.units['wave'])
    npix = int(np.round(np.log(wvmax/wvmin) / dlnlamb)) + 1
    new_wv = wvmin * np.exp(dlnlamb*np.arange(npix))
    # Final image
    f_flux = np.zeros((spec.nspec, npix))
    f_sig = np.zeros((spec.nspec, npix))
    f_wv = np.outer(np.ones(spec.nspec), new_wv.value)
    # Let's loop
    for ispec in range(spec.nspec):
        if debug:
            print("ispec={:d}".format(ispec))
        # Select
        spec.select = ispec
        # Rebin in obs frame
        tspec = spec.rebin(new_wv*(1+zarr[ispec]), do_sig=True, masking='none')
        # Save in rest-frame (worry about flambda)
        f_flux[ispec, :] = tspec.flux.value
        f_sig[ispec, :] = tspec.sig.value
    # Finish
    new_spec = XSpectrum1D(f_wv, f_flux, sig=f_sig, masking='none',
                           units=spec.units.copy())
    new_spec.meta = spec.meta.copy()
    # Return
    return new_spec


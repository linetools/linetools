""" Methods related to spectra
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import json
import warnings
import pdb

from astropy import units as u
from astropy import constants as const

from linetools import utils as liu


def meta_to_disk(in_meta):
    """ Polish up the meta dict for I/O

    Parameters
    ----------
    in_meta : dict
      list in 'headers' needs to be a series of Header objects
      or a series of dict objects

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
            if isinstance(header, dict):
                meta['headers'][kk] = header.copy()
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
      And to truncate original spectrum at
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
    gdp1 = np.where(spec1.wavelength <= wvmx)[0]
    gdp2 = np.where(spec2.wavelength > wvmx)[0]
    # Concatenate
    new_wv = np.concatenate((spec1.wavelength.value[gdp1],
                             spec2.wavelength.value[gdp2]))
    uwave = u.Quantity(new_wv, unit=spec1.units['wave'])
    new_fx = np.concatenate((spec1.flux.value[gdp1],
                             spec2.flux.value[gdp2] * scale))
    # Error
    if spec1.sig_is_set:
        new_sig = np.concatenate((spec1.sig[gdp1], spec2.sig[gdp2] * scale))
    else:
        new_sig = None

    # Continuum
    if spec1.co_is_set:
        new_co = np.concatenate((spec1.co[gdp1], spec2.co[gdp2] * scale))
    else:
        new_co = None

    # Generate
    spec3 = XSpectrum1D.from_tuple(
        (uwave, u.Quantity(new_fx), new_sig, new_co), meta=spec1.meta.copy())
    # Return
    return spec3


def collate(spectra):
    """ Generate a single XSpectrum1D instance containing an array of
    spectra from a list of individual XSpectrum1D spectra.
    Each spectrum is padded with extra pixels so that the
    wavelength ranges of all spectra are covered.
    Padded pixels are masked.

    Also note that masked pixels in the original data are ignored!

    Parameters
    ----------
    spectra : list
      of XSpectrum1D

    Returns
    -------
    new_spec : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    # Init
    maxpix = 0
    flg_co, flg_sig = False, False
    nspec = 0
    units = None
    # Setup
    for spec in spectra:
        nspec += spec.nspec
        for ii in range(spec.nspec):
            spec.select = ii
            maxpix = max(maxpix, spec.npix)
        if spec.co_is_set:
            flg_co = True
        if spec.sig_is_set:
            flg_sig = True
        # Check for identical units
        if units is None:
            units = spec.units
        else:
            assert spec.units == units
    # Generate data arrays
    wave = np.zeros((nspec, maxpix), dtype='float64')
    flux = np.zeros_like(wave, dtype='float32')
    if flg_sig:
        sig = np.zeros_like(wave, dtype='float32')
    else:
        sig = None
    if flg_co:
        co = np.zeros_like(flux)
    else:
        co = None
    # Load
    meta = dict(headers=[])
    idx = 0
    for xspec in spectra:
        # Allow for multiple spectra in the XSpectrum1D object
        for jj in range(xspec.nspec):
            xspec.select = jj
            wave[idx,:xspec.npix] = xspec.wavelength.value
            flux[idx,:xspec.npix] = xspec.flux.value
            if flg_sig:
                sig[idx,:xspec.npix] = xspec.sig.value
            if flg_co:
                if xspec.co_is_set:  # Allow for a mix of continua (mainly for specdb)
                    co[idx,:xspec.npix] = xspec.co.value
            idx += 1
        # Meta
        meta['headers'] += xspec.meta['headers']
    # Finish
    new_spec = XSpectrum1D(wave, flux, sig=sig, co=co, units=units.copy(),
                           masking='edges', meta=meta)
    # Return
    return new_spec


def rebin(spec, new_wv, do_sig=False, do_co=False, all=False, grow_bad_sig=False, **kwargs):
    """ Rebin a single spectrum in an XSpectrum1D object to a new wavelength array

    Uses simple linear interpolation.  The default (and only)
    option conserves counts (and flambda).

    WARNING: Do not trust either edge pixel of the new array.
      In fact the sig is set to 0 for each of these
    Also be aware that neighboring pixels are likely to be
    correlated in a manner that is not described by the error
    array.

    Parameters
    ----------
    new_wv : Quantity array
      New wavelength array
    do_sig : bool, optional
      Rebin error too (if it exists).
      S/N is only crudely conserved.
      Rejected pixels are propagated.
    do_co : bool, optional
      Rebin continuum if present
    all : bool, optional
      Rebin all spectra in the XSpectrum1D object?
    grow_bad_sig : bool, optional
      Allow sig<=0. values and grow them

    Returns
    -------
    newspec : XSpectrum1D
      XSpectrum1D of the rebinned spectrum
    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    from scipy.interpolate import interp1d
    # Save flux info to avoid unit issues
    funit = spec.flux.unit
    flux = spec.flux.value

    # Deal with nan
    badf = np.isnan(flux)
    if np.sum(badf) > 0:
        warnings.warn("Ignoring pixels with NAN in flux")
    gdf = ~badf
    flux = flux[gdf]

    # Check for bad pixels (not prepared for these)
    if spec.sig_is_set:
        sig = spec.sig.value
        bad_sig = sig[gdf] <= 0.
        if np.sum(bad_sig) > 0:
            if not grow_bad_sig:
                raise IOError("Data contains rejected pixels (sig=0). Use grow_bad_sig to proceed and grow them.")

    # Endpoints of original pixels
    npix = len(spec.wavelength)
    wvh = (spec.wavelength + np.roll(spec.wavelength, -1)) / 2.
    wvh[npix - 1] = spec.wavelength[npix - 1] + \
                    (spec.wavelength[npix - 1] - spec.wavelength[npix - 2]) / 2.
    dwv = wvh - np.roll(wvh, 1)
    dwv[0] = 2 * (wvh[0] - spec.wavelength[0])
    med_dwv = np.median(dwv.value)

    wvh = wvh[gdf]
    dwv = dwv[gdf]

    # Error
    if do_sig:
        if not spec.sig_is_set:
            raise IOError("sig must be set to rebin sig")
        var = sig[gdf]**2
        var[bad_sig] = 0.
    else:
        var = np.ones_like(flux)

    # Cumulative Sum
    cumsum = np.cumsum(flux * dwv)
    cumvar = np.cumsum(var * dwv)

    # Interpolate (loses the units)
    fcum = interp1d(wvh, cumsum, fill_value=0., bounds_error=False)
    fvar = interp1d(wvh, cumvar, fill_value=0., bounds_error=False)

    # Endpoints of new pixels
    nnew = len(new_wv)
    nwvh = (new_wv + np.roll(new_wv, -1)) / 2.
    nwvh[nnew - 1] = new_wv[nnew - 1] + \
                     (new_wv[nnew - 1] - new_wv[nnew - 2]) / 2.
    # Pad starting point
    bwv = np.zeros(nnew + 1) * new_wv.unit
    bwv[0] = new_wv[0] - (new_wv[1] - new_wv[0]) / 2.
    bwv[1:] = nwvh

    # Evaluate and put unit back
    newcum = fcum(bwv) * dwv.unit
    newvar = fvar(bwv) * dwv.unit

    # Endpoint
    #if (bwv[-1] > wvh[-1]):
    #    newcum[-1] = cumsum[-1]
    #    newvar[-1] = cumvar[-1]

    # Rebinned flux, var, co
    new_fx = (np.roll(newcum, -1) - newcum)[:-1]
    new_var = (np.roll(newvar, -1) - newvar)[:-1]

    # Normalize (preserve counts and flambda)
    new_dwv = bwv - np.roll(bwv, 1)
    #import pdb
    # pdb.set_trace()
    new_fx = new_fx / new_dwv[1:]
    # Preserve S/N (crudely)
    med_newdwv = np.median(new_dwv.value)
    new_var = new_var / (med_newdwv/med_dwv) / new_dwv[1:]

    # Return new spectrum
    if do_sig:
        # Create new_sig
        new_sig = np.zeros_like(new_var)
        gd = new_var > 0.
        new_sig[gd] = np.sqrt(new_var[gd].value)
        # Deal with bad pixels (grow_bad_sig should be True)
        bad = np.where(var <= 0.)[0]
        for ibad in bad:
            bad_new = np.where(np.abs(new_wv-spec.wavelength[ibad]) <
                               (new_dwv[1:]+dwv[ibad])/2)[0]
            new_sig[bad_new] = 0.
        # Zero out edge pixels -- not to be trusted
        igd = np.where(gd)[0]
        new_sig[igd[0]] = 0.
        new_sig[igd[-1]] = 0.
    else:
        new_sig = None

    # Continuum?
    if do_co:
        if not spec.co_is_set:
            raise IOError("Continuum must be set to request rebinning")
        co = spec.co.value
        co = co[gdf]
        cumco = np.cumsum(co * dwv)
        fco = interp1d(wvh, cumco, fill_value=0., bounds_error=False)
        newco = fco(bwv) * dwv.unit
        new_co = (np.roll(newco, -1) - newco)[:-1]
        new_co = new_co / new_dwv[1:]
    else:
        new_co = None

    # Finish
    newspec = XSpectrum1D.from_tuple((new_wv, new_fx*funit,
                                      new_sig, new_co),
                                     meta=spec.meta.copy(), **kwargs)
    # Return
    return newspec


def rebin_to_rest(spec, zarr, dv, debug=False, **kwargs):
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
    **kwargs :
      Passed to spec.rebin()

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
        tspec = spec.rebin(new_wv*(1+zarr[ispec]), do_sig=True, masking='none', **kwargs)
        # Save in rest-frame (worry about flambda)
        f_flux[ispec, :] = tspec.flux.value
        f_sig[ispec, :] = tspec.sig.value
    # Finish
    new_spec = XSpectrum1D(f_wv, f_flux, sig=f_sig, masking='none',
                           units=spec.units.copy())
    new_spec.meta = spec.meta.copy()
    # Return
    return new_spec


def smash_spectra(spec, method='average', debug=False):
    """ Collapse the data in XSpectrum1D.data
    One might call this 'stacking'
    Note: This works on the unmasked data array and returns
    unmasked spectra

    Parameters
    ----------
    spec : XSpectrum1D
    method : str, optional
      Approach to the smash ['smash']
    debug : bool, optional

    Returns
    -------
    new_spec : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    # Checks
    if spec.nspec <= 1:
        raise IOError("This method smashes an XSpectrum1D instance with multiple spectra")
    np.testing.assert_allclose(spec.data['wave'][0],spec.data['wave'][1])
    # Generate mask
    stack_msk = spec.data['sig'] > 0.
    if method == 'average':
        tot_flx = np.sum(spec.data['flux']*stack_msk,0)
        navg = np.sum(stack_msk,0)
        fin_flx = tot_flx / np.maximum(navg, 1.)
    else:
        raise IOError("Not prepared for this smash method")
    # Finish
    new_spec = XSpectrum1D(spec.data['wave'][0], fin_flx, masking='none',
                           units=spec.units.copy())
    new_spec.meta = spec.meta.copy()
    # Return
    return new_spec


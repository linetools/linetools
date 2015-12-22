""" Reading and writing of spectra
"""

from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str


# Import libraries
import numpy as np
import warnings
import os, pdb
import json

from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from astropy.table import Table, Column
from astropy.io.fits.hdu.table import BinTableHDU

#from xastropy.xutils import xdebug as xdb

def readspec(specfil, inflg=None, efil=None, verbose=False, flux_tags=None,
    sig_tags=None, multi_ivar=False, format='ascii', exten=None):
    """ Read a FITS file (or astropy Table or ASCII file) into a
    Spectrum1D class

    Parameters
    ----------
    specfil : str or Table
      Input file. If str:
        * FITS file must include in '.fit'
        * ASCII must either have a proper Table format
          or be 3 columns with WAVE,FLUX,ERROR
    efil : string, optional
      Explicit filename for Error array.  Code will attempt to find this
      file on its own.
    flux_tags : list of strings, optional
      Tags for flux in Binary FITS table. Default: flux_tags = ['SPEC',
      'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
    sig_tags : list of strings, optional
      Tags for error in Binary FITS table. Default : sig_tags =
      ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG', 'SIGMA_UP', 'ERRSTIS',
      'FLUXERR', 'er']
    multi_ivar : bool, optional
      If True, assume BOSS format of flux, ivar, log10(wave) in
      multi-extension FITS.
    format : str, optional
      Format for ASCII table input ['ascii']
    exten : int, optional
      FITS extension (mainly for multiple binary FITS tables)

    Returns
    -------
    An XSpectrum1D class
    """
    from specutils.io import read_fits as spec_read_fits
    from linetools.spectra.xspectrum1d import XSpectrum1D

    # Initialize

    if inflg == None:
        inflg = 0

    # Check specfil type
    if isinstance(specfil, Table):  # MAYBE SHOULD USE SPECUTILS FROM_TABLE
        datfil = 'None'
        # Dummy hdulist
        hdulist = [fits.PrimaryHDU(), specfil]
    elif isinstance(specfil, basestring):
        datfil = specfil.strip()
        flg_fits = False
        for ext in ['.fit']:
            if ext in specfil:
                flg_fits = True
        if flg_fits: # FITS
            # Read header
            datfil, chk = chk_for_gz(specfil.strip())
            if chk == 0:
                raise IOError('File does not exist {}'.format(specfil))
            hdulist = fits.open(os.path.expanduser(datfil))
        else: #ASCII
            tbl = Table.read(specfil,format=format)
            # No header?
            if tbl.colnames[0] == 'col1':
                names = 'WAVE', 'FLUX', 'ERROR', 'CONTINUUM'
                for i,name in enumerate(tbl.colnames):
                    tbl[name].name = names[i]
            hdulist = [fits.PrimaryHDU(), tbl]
    else:
        raise IOError('readspec: Bad spectra input')

    head0 = hdulist[0].header

    co = None

    ## #################
    # Binary FITS table?

    if is_UVES_popler(head0):
        xspec1d = parse_UVES_popler(hdulist)
    elif head0['NAXIS'] == 0:
        # Flux
        if flux_tags is None:
            flux_tags = ['SPEC', 'FLUX', 'FLAM', 'FX',
                         'FLUXSTIS', 'FLUX_OPT', 'fl', 'flux', 'counts']
        fx, fx_tag = get_table_column(flux_tags, hdulist, idx=exten)
        if fx is None:
            print('Binary FITS Table but no Flux tag')
            return
        # Error
        if sig_tags is None:
            sig_tags = ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG', 'SIGMA_UP',
                        'ERRSTIS', 'FLUXERR', 'sigma', 'sigma_flux', 'er','err']
        sig, sig_tag = get_table_column(sig_tags, hdulist)
        if sig is None:
            ivar_tags = ['IVAR', 'IVAR_OPT']
            ivar, ivar_tag = get_table_column(ivar_tags, hdulist, idx=exten)
            if ivar is None:
                var_tags = ['VAR', 'var']
                var, var_tag = get_table_column(var_tags, hdulist, idx=exten)
                sig = np.sqrt(var)
            else:
                sig = np.zeros(ivar.size)
                gdi = np.where( ivar > 0.)[0]
                sig[gdi] = np.sqrt(1./ivar[gdi])
        # Wavelength
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA','LOGLAM',
                     'WAVESTIS', 'WAVE_OPT', 'wa', 'wave']
        wave, wave_tag = get_table_column(wave_tags, hdulist, idx=exten)
        if wave_tag == 'LOGLAM':
            wave = 10.**wave
        if wave is None:
            print('Binary FITS Table but no wavelength tag')
            return
        co_tags = ['CONT', 'CO', 'CONTINUUM', 'co', 'cont']
        co, co_tag = get_table_column(co_tags, hdulist, idx=exten)

    elif head0['NAXIS'] == 1: # Data in the zero extension

        # How many entries?
        if len(hdulist) == 1: # Old school (one file per flux, error)
            # Error
            if efil == None:
                ipos = max(specfil.find('F.fits'),
                    specfil.find('f.fits'), specfil.find('flx.fits'))
                if ipos < 0:
                    # Becker XShooter style
                    ipos = specfil.find('.fits')
                    efil,chk = chk_for_gz(specfil[0:ipos]+'e.fits')
                else:
                    if specfil.find('F.fits') > 0:
                        efil,chk = chk_for_gz(specfil[0:ipos]+'E.fits')
                    else:
                        efil,chk = chk_for_gz(specfil[0:ipos]+'e.fits')
                    if efil is None:
                        efil,chk = chk_for_gz(specfil[0:ipos]+'err.fits')
                if efil is not None:
                    efil = os.path.expanduser(efil)

            # Error file
            if efil is not None:
                sig = fits.getdata(efil)
                uncertainty = StdDevUncertainty(sig)
            else:
                uncertainty = None

            #Log-Linear?
            try:
                dc_flag = head0['DC-FLAG']
            except KeyError:
                # The following is necessary for Becker's XShooter output
                cdelt1, dc_flag = get_cdelt_dcflag(head0)

            # Read
            if dc_flag == 0:
                # Read FITS file
                spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil), dispersion_unit='AA')
                spec1d.uncertainty = uncertainty
                xspec1d = XSpectrum1D.from_spec1d(spec1d)
            elif dc_flag == 1: # Generate wavelengths and use array approach
                fx = hdulist[0].data
                # Generate wave
                wave = setwave(head0)
            else:
                raise ValueError('DC-FLAG has unusual value {:d}'.format(dc_flag))

        elif hdulist[0].name == 'FLUX':
            # NEW SCHOOL (one file for flux and error)
            if 'WAVELENGTH' not in hdulist:
                spec1d = spec_read_fits.read_fits_spectrum1d(
                    os.path.expanduser(datfil), dispersion_unit='AA')
                xspec1d = XSpectrum1D.from_spec1d(spec1d)
            else:
                wave = hdulist['WAVELENGTH'].data * u.AA
                fx = hdulist['FLUX'].data
                xspec1d = XSpectrum1D.from_array(wave, u.Quantity(fx))

            # Error array
            if 'ERROR' in hdulist:
                sig = hdulist['ERROR'].data
                xspec1d.uncertainty = StdDevUncertainty(sig)
            else:
                sig = None

            if 'CONTINUUM' in hdulist:
                xspec1d.co = hdulist['CONTINUUM'].data

            if 'METADATA' in head0:
                xspec1d.meta.update(json.loads(head0['METADATA']))

        else:  # ASSUMING MULTI-EXTENSION
            if len(hdulist) <= 2:
                raise RuntimeError('No wavelength info but only 2 extensions!')
            fx = hdulist[0].data.flatten()
            try:
                sig = hdulist[1].data.flatten()
            except AttributeError:  # Error array is "None"
                sig = None
            wave = hdulist[2].data.flatten()
            # BOSS/SDSS?
            try:
                multi_ivar = head0['TELESCOP'][0:4] in ['SDSS']
            except KeyError:
                pass
            #
            if multi_ivar is True:
                tmpsig = np.zeros(len(sig))
                gdp = np.where(sig > 0.)[0]
                tmpsig[gdp] = np.sqrt(1./sig[gdp])
                sig = tmpsig
                wave = 10.**wave
    elif (head0['NAXIS'] == 2) and (head0['NAXIS2'] == 5): # SDSS .fit format
        fx = hdulist[0].data[0,:].flatten()
        sig = hdulist[0].data[2,:].flatten()
        wave = setwave(head0)
    else:  # Should not be here
        print('Looks like an image')
        return


    # Generate, as needed
    if 'xspec1d' not in locals():
        # Give Ang as default
        if not hasattr(wave, 'unit'):
            uwave = u.Quantity(wave, unit=u.AA)
        elif wave.unit is None:
            uwave = u.Quantity(wave, unit=u.AA)
        else:
            uwave = u.Quantity(wave)
        if sig is not None:
            xspec1d = XSpectrum1D.from_array(uwave, u.Quantity(fx),
                                             uncertainty=StdDevUncertainty(sig))
        else:
            xspec1d = XSpectrum1D.from_array(uwave, u.Quantity(fx))

    if np.any(np.isnan(xspec1d.dispersion)):
        warnings.warn('WARNING: Some wavelengths are NaN')

    # Filename
    xspec1d.filename = datfil

    if not hasattr(xspec1d, 'co'):
        xspec1d.co = co
        # Final check for continuum in a separate file
        if isinstance(specfil,basestring):
            if co is None and specfil.endswith('.fits'):
                try:
                    xspec1d.co = fits.getdata(specfil.replace('.fits', '_c.fits'))
                except IOError:
                    pass

    # Add in the header
    xspec1d.head = head0

    # Return
    return xspec1d


#### ###############################
#  Grab values from the Binary FITS Table or Table
def get_table_column(tags, hdulist, idx=None):
    """ Find a column in a FITS binary table

    Used to return flux/error/wave values from a binary FITS table
    from a list of tags.

    Parameters
    ----------
    tags : list
     List of string tag names
    hdulist : fits header data unit list  
    idx : int, optional
     Index of list for Table input

    Returns
    -------
    dat : float array
      Data values corresponding to the first tag found
      Returns None if no match
    """
    if idx is None:
        idx = 1
    dat = None
    # Use Table
    if isinstance(hdulist[idx],BinTableHDU):
        tab = Table(hdulist[idx].data)
    else:
        tab = hdulist[idx]

    # Grab
    names = set(tab.dtype.names)
    for tag in tags:
        if tag in names:
            dat = tab[tag]
            break  # Break with first hit

    # Return
    if dat is not None:
        return dat.flatten(), tag
    else:
        return dat, 'NONE'

#### ###############################
#  Set wavelength array using Header cards
def setwave(hdr):
    ''' Generate wavelength array from a header

    Parameters
    ----------
    hdr : FITS header

    Returns
    -------
    wave : ndarray
      No units yet
    '''

    # Parse the header
    npix = hdr['NAXIS1']
    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else 1.
    crval1 = hdr['CRVAL1']

    cdelt1, dc_flag = get_cdelt_dcflag(hdr)

    # Generate
    wave = crval1 + cdelt1 * (np.arange(npix) + 1. - crpix1)
    if dc_flag == 1:
        wave = 10.**wave # Log

    return wave

def get_cdelt_dcflag(hd):
    """ Find the wavelength stepsize and dcflag from a fits header.

    Parameters
    ----------
    hd : astropy.io.fits header instance

    Returns
    -------
    cdelt, dc_flag : float, int
      Wavelength stepsize and dcflag (1 if log-linear scale, 0 if linear).
    """
    cdelt = None
    if 'CDELT1' in hd:
        cdelt1 = hd['CDELT1']
    elif 'CD1_1' in hd:
        cdelt1 = hd['CD1_1']  # SDSS style

    dc_flag = 0
    if 'DC-FLAG' in hd:
        dc_flag = hd['DC-FLAG']
    elif cdelt1 < 1e-4:
        import warnings
        warnings.warn('WARNING: CDELT1 < 1e-4, Assuming log wavelength scale')
        dc_flag = 1

    return cdelt1, dc_flag


#### ###############################
# Deal with .gz extensions, usually on FITS files
#   See if filenm exists, if so pass it back
#
def chk_for_gz(filenm):
    """ Checks for .gz extension to an input filename and returns file

    Also parses the ~ if given

    Parameters
    ----------
    filenm : string
     Filename to query

    Returns
    -------
    filenm+XX : string
      Returns in this order:
        i. Input filename if it exists
        ii. Input filename if it has .gz extension already
        iii. Input filename.gz if that exists
        iv. Input filename.gz if that exists
    chk : bool or int
      * True if file exists
      * 0 if No check was performed
      * False if no file exists
    """
    import os
    from os.path import expanduser
    filenm = expanduser(filenm)

    # File exist?
    if os.path.lexists(filenm):
        chk=True
        return filenm, chk

    # .gz already
    if filenm.find('.gz') > 0:
        chk=0
        return filenm, chk

    # Add .gz
    if os.path.lexists(filenm+'.gz'):
        chk=True
        return filenm+'.gz', chk
    else:
        chk=False
        return None, chk

def is_UVES_popler(hd):
    """ Check if this header is UVES_popler output.
    """
    if 'history' not in hd:
        return False
    for row in hd['history']:
        if 'UVES POst Pipeline Echelle Reduction' in row:
            return True
    return False

def parse_UVES_popler(hdulist):
    """ Read a spectrum from a UVES_popler-style fits file.
    """
    from linetools.spectra.xspectrum1d import XSpectrum1D

    hd = hdulist[0].header
    uwave = setwave(hd) * u.Angstrom
    co = hdulist[0].data[3]
    fx = hdulist[0].data[0] * co  #  Flux
    sig = hdulist[0].data[1] * co
    xspec1d = XSpectrum1D.from_array(uwave, u.Quantity(fx),
                                     uncertainty=StdDevUncertainty(sig))
    xspec1d.co = co
    return xspec1d

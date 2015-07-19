"""
Module for read/write of spectra FITS files
  -- Offers greater flexibility than the code in specutils
"""

from __future__ import print_function, absolute_import, division, unicode_literals

# Import libraries
import numpy as np
import os

from astropy.io import fits, ascii
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from astropy.table import Table, Column
from astropy.io.fits.fitsrec import FITS_rec
from astropy.io.fits.hdu.table import BinTableHDU

#from xastropy.xutils import xdebug as xdb
from linetools.spectra.utils import XSpectrum1D

#### ###############################
#  Generate Spectrum1D from FITS file
#
def readspec(specfil, inflg=None, efil=None, verbose=False, flux_tags=None, 
    sig_tags=None, multi_ivar=False):
    ''' Read a FITS file (or astropy Table or ASCII file) into a Spectrum1D class

    Parameters:
    -----------
    specfil: str or Table
      Input file
      If str, 
        FITS file must include in '.fit'
        ASCII must either have a proper Table format
          or be 3 columns with WAVE,FLUX,ERROR
    efil: string, optional
      Explicit filename for Error array.  Code will attempt to find this
      file on its own.
    flux_tags: list of strings, optional
      Tags for flux in Binary FITS table
      Default: flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
    sig_tags: list of strings, optional
      Tags for error in Binary FITS table
      Default: flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
    multi_ivar: bool, optional 
      If True, assume BOSS format of  flux, ivar, log10(wave) in multi-extension FITS

    Returns:
    -----------
    A Spectrum1D class (or XSpectrum1D which is an overloaded variant)   
    '''
    from specutils.io import read_fits as spec_read_fits

    # Initialize
    dat = None
    if inflg == None:
        inflg = 0

    # Check specfil type
    if isinstance(specfil,Table):
        datfil = 'None'
        # Dummy hdulist
        hdulist = [fits.PrimaryHDU(), specfil]
    elif isinstance(specfil, basestring): 
        flg_fits = False
        for ext in ['.fit']:
            if ext in specfil:
                flg_fits = True
        if flg_fits: # FITS
            # Read header
            datfil,chk = chk_for_gz(specfil)
            if chk == 0:
                print('xastropy.spec.readwrite: File does not exist ', specfil)
                return -1
            hdulist = fits.open(os.path.expanduser(datfil))
        else: #ASCII
            try:
                tbl = Table.read(specfil)
            except Exception:
                tbl = ascii.read(specfil, names=['WAVE', 'FLUX', 'ERROR'])
            hdulist = [fits.PrimaryHDU(), tbl]
    else:
        raise IOError('readspec: Bad spectra input')

    head0 = hdulist[0].header

    ## #################
    # Binary FITS table?
    if head0['NAXIS'] == 0:
        # Flux 
        if flux_tags is None:
            flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
        fx, fx_tag = get_table_column(flux_tags, hdulist)
        if fx is None:
            print('spec.readwrite: Binary FITS Table but no Flux tag')
            return
        # Error
        if sig_tags is None:
            sig_tags = ['ERROR','ERR','SIGMA_FLUX','FLAM_SIG', 'SIGMA_UP', 'ERRSTIS', 'FLUXERR', 'er']
        sig, sig_tag = get_table_column(sig_tags, hdulist)
        if sig is None:
            ivar_tags = ['IVAR', 'IVAR_OPT']
            ivar, ivar_tag = get_table_column(ivar_tags, hdulist)
            if ivar is None:
                print('spec.readwrite: Binary FITS Table but no error tags')
                return
            else: 
                sig = np.zeros(ivar.size)
                gdi = np.where( ivar > 0.)[0]
                sig[gdi] = np.sqrt(1./ivar[gdi])
        # Wavelength
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA','LOGLAM', 'WAVESTIS', 'WAVE_OPT', 'wa']
        wave, wave_tag = get_table_column(wave_tags, hdulist)
        if wave_tag == 'LOGLAM':
            wave = 10.**wave
        if wave is None:
            print('spec.readwrite: Binary FITS Table but no wavelength tag')
            return
    elif head0['NAXIS'] == 1: # Data in the zero extension
        # How many entries?
        if len(hdulist) == 1: # Old school (one file per flux, error)
            # Error
            if efil == None:
                ipos = max(specfil.find('F.fits'),specfil.find('f.fits'))
                if ipos < 0: # No error array
                    efil = None
                    #sig = np.zeros(fx.size)
                else:
                    if specfil.find('F.fits') > 0:
                        efil,chk = chk_for_gz(specfil[0:ipos]+'E.fits')
                    else:
                        efil,chk = chk_for_gz(specfil[0:ipos]+'e.fits')
                if efil != None:
                    efil=os.path.expanduser(efil)

            # Error file
            if efil != None:
                # BZERO -- This has to happen before reading the data!
                try:
                    bzero = head0['BZERO']
                except KeyError:
                    bzero = 0.
                sig=fits.getdata(efil) - bzero
                uncertainty = StdDevUncertainty(sig)
            else:
                uncertainty = None

            #Log-Linear?
            try:
                dc_flag = head0['DC-FLAG']
            except KeyError:
                dc_flag = 0
 
            if dc_flag == 0:
                # Read FITS file
                spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil),
                                                             dispersion_unit='AA')
                spec1d.uncertainty = uncertainty
                xspec1d = XSpectrum1D.from_spec1d(spec1d) # UNTESTED!
            elif dc_flag == 1: # Generate wavelengths and use array approach
                fx = hdulist[0].data
                # Generate wave
                wave = setwave(head0)
            else:
                raise ValueError('DC-FLAG has unusual value {:d}'.format(dc_flag))


        elif len(hdulist) == 2: # NEW SCHOOL (one file per flux, error)
            spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil), dispersion_unit='AA')
            # Error array
            sig = hdulist[1].data
            spec1d.uncertainty = StdDevUncertainty(sig)
            #
            xspec1d = XSpectrum1D.from_spec1d(spec1d)

        else:  # ASSUMING MULTI-EXTENSION
            if len(hdulist) <= 2:
                print('spec.readwrite: No wavelength info but only 2 extensions!')
                return
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
    else:  # Should not be here
        print('spec.readwrite: Looks like an image')
        return dat

    # Generate, as needed
    if 'xspec1d' not in locals():
        # Give Ang as default
        if not hasattr(wave, 'unit'):
            uwave = u.Quantity(wave, unit=u.AA)
        else:
            if wave.unit is None:
                uwave = u.Quantity(wave, unit=u.AA)
            else:
                uwave = u.Quantity(wave)
        xspec1d = XSpectrum1D.from_array(uwave, u.Quantity(fx),
                                         uncertainty=StdDevUncertainty(sig))

    xspec1d.filename = specfil

    # Continuum?
    try:
        co = fits.getdata(name+'_c.fits')
    except:
        try:
            npix = len(fx)
        except UnboundLocalError:
            npix = len(xspec1d.flux)
        co = np.nan*np.ones(npix)

    # Add in the header
    xspec1d.head = head0

    # Return 
    return xspec1d


#### ###############################
#  Grab values from the Binary FITS Table or Table
def get_table_column(tags, hdulist, idx=1):
    '''Simple script to return values from a FITS table
    Function used to return flux/error/wave values from 
    a binary FITS table from a list of tags

    Parameters:
    -----------
    tags: list
     List of string tag names
    idx: int, optional
     Index of list for Table input [Default: 1]

    Returns:
    -----------
    dat: float array
      Data values corresponding to the first tag found
      Returns None if no match
    '''
    dat = None
    # Use Table
    if isinstance(hdulist[idx],BinTableHDU):
        tab = Table(hdulist[idx].data)
    else:
        tab = hdulist[idx]

    # Grab
    for tag in tags:
        if tag in tab.dtype.names: 
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
    ''' Generate wavelength array from header
    '''

    # Parse the header
    npix = hdr['NAXIS1'] 
    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else 1.
    crval1 = hdr['CRVAL1']
    cdelt1 = hdr['CDELT1']
    dcflag = hdr['DC-FLAG'] if 'DC-FLAG' in hdr else None

    if cdelt1 < 1e-4:
        import warnings
        warnings.warn('WARNING: CDELT1 < 1e-4, Assuming log wavelength scale')
        dcflag = 1

    # Generate
    wave = crval1 + cdelt1 * (np.arange(npix) + 1. - crpix1)
    if dcflag == 1:
        wave = 10.**wave # Log

    return wave

#### ###############################
# Deal with .gz extensions, usually on FITS files
#   See if filenm exists, if so pass it back
#
def chk_for_gz(filenm):
    '''Checks for .gz extension to an input filename and returns file
    Also parses the ~ if given

    Parameters:
    -----------
    filenm: string
     Filename to query

    Returns:
    -----------
    filenm+XX: string
      Returns in this order:
        i. Input filename if it exists
        ii. Input filename if it has .gz extension already
        iii. Input filename.gz if that exists
        iv. Input filename.gz if that exists
    chk: bool or int
      True if file exists
      0 if No check was performed
      False if no file exists 
    '''
    import os, pdb
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

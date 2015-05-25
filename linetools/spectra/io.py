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

from specutils.io import read_fits as spec_read_fits

from xastropy.xutils import xdebug as xdb
from linetools.spectra.utils import XSpectrum1D

#### ###############################
#  Generate Spectrum1D from FITS file
#
def readspec(specfil, inflg=None, efil=None, verbose=False, flux_tags=None, 
    sig_tags=None, multi_ivar=False):
    ''' Read a FITS file (or astropy Table) into a Spectrum1D class

    Parameters:
    -----------
    specfil: string or Table
      Input file
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

    # Initialize
    dat = None
    if inflg == None:
        inflg = 0

    # Check specfil type
    if type(specfil) is Table:
        datfil = 'None'
        # Dummy hdulist
        hdulist = [fits.PrimaryHDU(), specfil]
    else:
        # Read header
        datfil,chk = chk_for_gz(specfil)
        if chk == 0:
            print('xastropy.spec.readwrite: File does not exist ', specfil)
            return -1
        hdulist = fits.open(os.path.expanduser(datfil))

    head0 = hdulist[0].header

    ## #################
    # Binary FITS table?
    if head0['NAXIS'] == 0:
        # Flux 
        if flux_tags is None:
            flux_tags = ['SPEC', 'FLUX','FLAM','FX', 'FLUXSTIS', 'FLUX_OPT', 'fl']
        fx, fx_tag = get_table_column(flux_tags, hdulist)
        #xdb.set_trace()
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
            # Generate Spectrum1D
            spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil),
                                                         dispersion_unit='AA',
                                                         efil=efil)
            xspec1d = XSpectrum1D.from_spec1d(spec1d)

            #spec1d = spec_read_fits.read_fits_spectrum1d(os.path.expanduser(datfil))

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
            sig = hdulist[1].data.flatten()
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
#### ###############################
#  Set wavelength array using Header cards
def setwave(hdr):

    # DEPRECATED
    xdb.set_trace()
    # Initialize
    SCL = 1.
    
    # Parse the header
    npix = hdr['NAXIS1'] 
    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else 1.
    crval1 = hdr['CRVAL1'] if 'CRVAL1' in hdr else 1.
    cdelt1 = hdr['CDELT1'] if 'CDELT1' in hdr else 1.
    ctype1 = hdr['CTYPE1'] if 'CTYPE1' in hdr else None
    dcflag = hdr['DC-FLAG'] if 'DC-FLAG' in hdr else None

    # Generate
    if (dcflag == 1) or (cdelt1 < 1e-4):
        wave = SCL * 10.**(crval1 + ( cdelt1 * np.arange(npix) + 1. - crpix1) ) # Log
    xdb.set_trace()

    # Return
    return wave

#### ###############################
#### ###############################
#  Grab values from the Binary FITS Table or Table
def get_table_column(tags, hdulist):
    dat = None
    ii = 0
    # Use Table
    if type(hdulist[1]) is BinTableHDU:
        tab = Table(hdulist[1].data)
    else:
        tab = hdulist[1]

    # Grab
    for tag in tags:
        if tag in tab.dtype.names: 
            dat = tab[tag]
            break  # Break with first hit

    '''
    For BinTableHDU (deprecated)
    while(ii < len(tags)):
        if tags[ii] in hdulist[1].columns.names: 
            dat = hdulist[1].data[tags[ii]]
            break  # Break with first hit
        else:
            ii = ii + 1
    '''
    # Return
    if dat is not None:
        return dat.flatten(), tag
    else: 
        return dat, 'NONE'


#### ###############################
# Deal with .gz extensions, usually on FITS files
#   See if filenm exists, if so pass it back
#
def chk_for_gz(filenm,chk=None):

    import os, pdb

    # File exist?
    if os.path.lexists(filenm): 
        chk=1
        return filenm, chk

    # .gz already
    if filenm.find('.gz') > 0:
        chk=0
        return filenm, chk

    # Add .gz
    if os.path.lexists(filenm+'.gz'): 
        chk=1
        return filenm+'.gz', chk
    else:
        chk=0
        return filenm, chk







#### ###############################
# Testing
if __name__ == '__main__':
    flg_test = 0
    flg_test += 1 # MagE
    flg_test += 2**1 # LRIS LowRedux

    # Standard log-linear read (MagE)
    if (flg_test % 2**1) >= 2**0:
        fil = '~/PROGETTI/LLSZ3/data/normalize/UM669_nF.fits'
        #fil = '/Users/xavier/Dropbox/QSOPairs/data/MAGE_redux/SDSSJ085357.49-001106.1_F.fits.gz'
        #efil = '~ers/xavier/PROGETTI/LLSZ3/data/normalize/UM669_nE.fits'
        myspec = readspec(fil)
        #xdb.xplot(myspec.dispersion, myspec.flux)
        myspec.plot()
        #xdb.xplot(myspec.dispersion, myspec.flux, myspec.uncertainty.array)

    # LowRedux
    if (flg_test % 2**2) >= 2**1:
        fil = '/Users/xavier/Dropbox/QSOPairs/data/LRIS_redux/SDSSJ231254.65-025403.1_b400_F.fits.gz'
        myspec = readspec(fil)
        xdb.xplot(myspec.dispersion, myspec.flux, myspec.uncertainty.array)
        #xdb.set_trace()

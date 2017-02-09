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
from astropy import units as u
from astropy.table import Table, Column
from astropy.io.fits.hdu.table import BinTableHDU

from .xspectrum1d import XSpectrum1D


def readspec(specfil, inflg=None, efil=None, verbose=False, multi_ivar=False,
             format='ascii', exten=None, head_exten=0, debug=False, select=0,
             **kwargs):
    """ Read a FITS file (or astropy Table or ASCII file) into a
    XSpectrum1D class

    Parameters
    ----------
    specfil : str or Table
      Input file. If str:
        * FITS file are detected by searching for '.fit' in their filename.
        * ASCII must either have a proper Table format or be 3 (WAVE,
          FLUX, ERROR) or 4 (WAVE, FLUX, ERROR, CONTINUUM) columns. If
          the file has more than 4 columns with no header it will raise an error.
    efil : string, optional
      A filename for Error array, if it's in a separate file to the
      flux. The code will attempt to find this file on its own.
    multi_ivar : bool, optional
      If True, assume BOSS format of flux, ivar, log10(wave) in a
      multi-extension FITS.
    format : str, optional
      Format for ASCII table input. Default 'ascii'.
    exten : int, optional
      FITS extension (mainly for multiple binary FITS tables)
    select : int, optional
      Selected spectrum (for sets of 1D spectra, e.g. DESI brick)
    head_exten : int, optional
      Extension for header to ingest

    Returns
    -------
    An XSpectrum1D class

    """

    # Initialize
    if inflg is None:
        inflg = 0

    # Check specfil type
    if isinstance(specfil, Table):
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
        elif '.hdf5' in specfil:  # HDF5
            return parse_hdf5(specfil, **kwargs)
        else: #ASCII
            tbl = Table.read(specfil,format=format)
            # No header?
            if tbl.colnames[0] == 'col1':
                if len(tbl.colnames) > 4:
                    raise IOError('No header found in ASCII file {}, \
                                  and has more than four columns. Please check its format.'.format(specfil))
                names = ['WAVE', 'FLUX', 'ERROR', 'CONTINUUM']
                for i, name in enumerate(tbl.colnames):
                    tbl[name].name = names[i]
                warnings.warn('No header found in ASCII file {}, assuming columns to be: {}'.format(specfil, names[:len(tbl.colnames)]))
            # import pdb; pdb.set_trace()
            hdulist = [fits.PrimaryHDU(), tbl]
    else:
        raise IOError('readspec: Bad spectra input.')

    head0 = hdulist[0].header

    if is_UVES_popler(head0):
        if debug:
            print('linetools.spectra.io.readspec(): Reading UVES popler format')
        xspec1d = parse_UVES_popler(hdulist)

    elif head0['NAXIS'] == 0:
        # Binary FITS table
        if debug:
            print('linetools.spectra.io.readspec(): Assuming binary fits table')
        xspec1d = parse_FITS_binary_table(hdulist, exten=exten, **kwargs)

    elif head0['NAXIS'] == 1: # Data in the zero extension

        # How many entries?
        if len(hdulist) == 1:  # Old school (one file per flux, error)
            if debug:
                print(
  'linetools.spectra.io.readspec(): Assuming flux and err in separate files')
            xspec1d = parse_two_file_format(specfil, hdulist, efil=efil)

        elif hdulist[0].name == 'FLUX':
            if debug:
                print(
  'linetools.spectra.io.readspec(): Assuming separate flux and err files.')
            xspec1d = parse_linetools_spectrum_format(hdulist)

        else:  # ASSUMING MULTI-EXTENSION
            co=None
            if debug:
                print(
              'linetools.spectra.io.readspec(): Assuming multi-extension')

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

            # Look for co
            if len(hdulist) == 4:
                data = hdulist[3].data
                if 'float' in data.dtype.name:  # This can be an int mask (e.g. BOSS)
                    co = data

            wave = give_wv_units(wave)
            xspec1d = XSpectrum1D.from_tuple((wave, fx, sig, co), **kwargs)

    elif head0['NAXIS'] == 2:
        if (hdulist[0].name == 'FLUX') and (hdulist[2].name == 'WAVELENGTH'):  # DESI
            if debug:
                print('linetools.spectra.io.readspec(): Assuming DESI brick')
            xspec1d = parse_DESI_brick(hdulist, select=select)

        else:  # SDSS
            if debug:
                print('linetools.spectra.io.readspec(): Assuming SDSS format')
            fx = hdulist[0].data[0, :].flatten()
            sig = hdulist[0].data[2, :].flatten()
            wave = setwave(head0)
            xspec1d = XSpectrum1D.from_tuple(
                (give_wv_units(wave), fx, sig, None))
    else:  # Should not be here
        print('Not sure what has been input.  Send to JXP.')
        return


    # Generate, as needed
    """
    if 'xspec1d' not in locals():
        # Give Ang as default
        if not hasattr(wave, 'unit'):
            uwave = u.Quantity(wave, unit=u.AA)
        elif wave.unit is None:
            uwave = u.Quantity(wave, unit=u.AA)
        else:
            uwave = u.Quantity(wave)
        xspec1d = XSpectrum1D.from_tuple((wave, fx, sig, None))
    """

    if np.any(np.isnan(xspec1d.wavelength)):
        warnings.warn('WARNING: Some wavelengths are NaN')

    # Filename
    xspec1d.filename = datfil

    if not xspec1d.co_is_set:
        # Final check for continuum in a separate file
        if isinstance(specfil, basestring) and (specfil.endswith('.fits') or specfil.endswith('.fits.gz')):
            co_filename = specfil.replace('.fits', '_c.fits')
            if os.path.exists(co_filename):
                tmpco = fits.getdata(co_filename)
                if tmpco.size != xspec1d.totpix:
                    warnings.warn("Continuum size does not match native spectrum")
                    warnings.warn("Continuing under the assumption that this is due to a masked array")
                    gdp = ~xspec1d.data['flux'][xspec1d.select].mask
                    xspec1d.data['co'][xspec1d.select][gdp] = tmpco
                else:
                    xspec1d.data['co'][xspec1d.select] = tmpco
                # Mask
                xspec1d.data['co'][xspec1d.select].mask = xspec1d.data['flux'][xspec1d.select].mask

    # Add in the header
    if head_exten == 0:
        xspec1d.meta['headers'][0] = head0
    else:
        head = hdulist[head_exten].header
        xspec1d.meta['headers'][0] = head
    if xspec1d.nspec > 1:
        warnings.warn("Read in only 1 header (into meta['headers'][0]")

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
            dat = np.array(tab[tag])
            break  # Break with first hit

    # Return
    if dat is not None:
        return dat.flatten(), tag
    else:
        return dat, 'NONE'


def get_wave_unit(tag, hdulist, idx=None):
    """ Attempt to pull wavelength unit from the Table
    Parameters
    ----------
    tag : str
     Tag used for wavelengths
    hdulist : fits header data unit list
    idx : int, optional
     Index of list for Table input

    Returns
    -------
    unit : astropy Unit
      Defaults to None
    """
    from astropy.units import Unit
    if idx is None:
        idx = 1
    # Use Table
    if isinstance(hdulist[idx],BinTableHDU):
        tab = Table(hdulist[idx].data)
        header = hdulist[idx].header
    else:
        # NEED HEADER INFO
        return None
    # Try table header (following VLT/X-Shooter here)
    keys = header.keys()
    values = header.values()
    hidx = values.index(tag)
    if keys[hidx][0:5] == 'TTYPE':
        try:
            tunit = header[keys[hidx].replace('TYPE','UNIT')]
        except KeyError:
           return None
        else:
            if tunit in ['Angstroem', 'Angstroms']:
                tunit = 'Angstrom'
            unit = Unit(tunit)
            return unit
    else:
        return None


#### ###############################
#  Set wavelength array using Header cards
def setwave(hdr):
    """ Generate wavelength array from a header

    Parameters
    ----------
    hdr : FITS header

    Returns
    -------
    wave : ndarray
      No units yet
    """

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


def give_wv_units(wave):
    """ Give a wavelength array units of Angstroms, if unitless.

    Parameters
    ----------
    wave : array or Quantity
      Input wavelength array

    Returns
    -------
    uwave: Quantity
      Output wavelengths in Angstroms if input is unitless, or the
      input array unchanged otherwise.
    """
    if not hasattr(wave, 'unit'):
        uwave = u.Quantity(wave, unit=u.AA)
    elif wave.unit is None:
        uwave = u.Quantity(wave, unit=u.AA)
    else:
        uwave = u.Quantity(wave)

    return uwave


def is_UVES_popler(hd):
    """ Check if this header is UVES_popler output.

    Parameters
    ----------
    hd : FITS header

    Returns
    -------
    True if a UVES_popler file, False otherwise.
    """
    if 'history' not in hd:
        return False
    for row in hd['history']:
        if 'UVES POst Pipeline Echelle Reduction' in row:
            return True
    return False

def parse_UVES_popler(hdulist):
    """ Read a spectrum from a UVES_popler-style fits file.

    Parameters
    ----------
    hdulist : FITS HDU list

    Returns
    -------
    xspec1d : XSpectrum1D
      Parsed spectrum
    """
    from linetools.spectra.xspectrum1d import XSpectrum1D

    hd = hdulist[0].header
    uwave = setwave(hd) * u.Angstrom
    co = hdulist[0].data[3]
    fx = hdulist[0].data[0] * co  #  Flux
    sig = hdulist[0].data[1] * co
    xspec1d = XSpectrum1D.from_tuple((uwave, fx, sig, co))
    return xspec1d

def parse_FITS_binary_table(hdulist, exten=None, wave_tag=None, flux_tag=None,
                            sig_tag=None, co_tag=None, var_tag=None,
                            ivar_tag=None, **kwargs):
    """ Read a spectrum from a FITS binary table

    Parameters
    ----------
    hdulist : FITS HDU list
    exten : int, optional
      Extension for the binary table.
    wave_tag : str, optional
    flux_tag : str, optional
    sig_tag : str, optional
    co_tag : str, optional
    var_tag : str, optional

    Returns
    -------
    xspec1d : XSpectrum1D
      Parsed spectrum
    """
    # Flux
    if flux_tag is None:
        flux_tags = ['SPEC', 'FLUX', 'FLAM', 'FX', 'FNORM',
                     'FLUXSTIS', 'FLUX_OPT', 'fl', 'flux', 'counts',
                     'COUNTS']
    else:
        flux_tags = [flux_tag]
    fx, fx_tag = get_table_column(flux_tags, hdulist, idx=exten)
    if fx is None:
        print('Binary FITS Table but no Flux tag. Searched fo these tags:\n',
              flux_tags)
        return
    # Error
    if sig_tag is None:
        sig_tags = ['ERROR','ERR','SIGMA_FLUX','ERR_FLUX', 'ENORM', 'FLAM_SIG', 'SIGMA_UP',
                    'ERRSTIS', 'FLUXERR', 'SIGMA', 'sigma', 'sigma_flux',
                    'er', 'err', 'error', 'sig', 'fluxerror']
    else:
        sig_tags = [sig_tag]
    sig, sig_tag = get_table_column(sig_tags, hdulist)
    if sig is None:
        if ivar_tag is None:
            ivar_tags = ['IVAR', 'IVAR_OPT', 'ivar', 'FLUX_IVAR']
        else:
            ivar_tags = [ivar_tag]
        ivar, ivar_tag = get_table_column(ivar_tags, hdulist, idx=exten)
        if ivar is None:
            if var_tag is None:
                var_tags = ['VAR', 'var']
            else:
                var_tags = [var_tag]
            var, var_tag = get_table_column(var_tags, hdulist, idx=exten)
            if var is None:
                warnings.warn('No error tag found. Searched for these tags:\n'+ str(sig_tags + ivar_tags + var_tags))
            else:
                sig = np.sqrt(var)
        else:
            sig = np.zeros(ivar.size)
            gdi = np.where( ivar > 0.)[0]
            sig[gdi] = np.sqrt(1./ivar[gdi])
    # Wavelength
    if wave_tag is None:
        wave_tags = ['WAVE','WAVELENGTH','LAMBDA','LOGLAM',
                     'WAVESTIS', 'WAVE_OPT', 'wa', 'wave', 'loglam','wl']
    else:
        wave_tags = [wave_tag]
    wave, wave_tag = get_table_column(wave_tags, hdulist, idx=exten)
    if wave_tag in ['LOGLAM','loglam']:
        wave = 10.**wave
    # Try for unit
    wv_unit = get_wave_unit(wave_tag, hdulist, idx=exten)
    if wv_unit is not None:
        wave = wave * wv_unit
    if wave is None:
        print('Binary FITS Table but no wavelength tag. Searched for these tags:\n',
              wave_tags)
        return
    if co_tag is None:
        co_tags = ['CONT', 'CO', 'CONTINUUM', 'co', 'cont']
    else:
        co_tags = [co_tag]
    co, co_tag = get_table_column(co_tags, hdulist, idx=exten)
    # Finish
    xspec1d = XSpectrum1D.from_tuple((give_wv_units(wave), fx, sig, co), **kwargs)

    if 'METADATA' in hdulist[0].header:
        xspec1d.meta.update(json.loads(hdulist[0].header['METADATA']))
    return xspec1d

def parse_linetools_spectrum_format(hdulist):
    """ Parse an old linetools-format spectrum from an hdulist

    Parameters
    ----------
    hdulist : FITS HDU list

    Returns
    -------
    xspec1d : XSpectrum1D
      Parsed spectrum

    """
    if 'WAVELENGTH' not in hdulist:
        pdb.set_trace()
        xspec1d = XSpectrum1D.from_spec1d(spec1d)
    else:
        wave = hdulist['WAVELENGTH'].data * u.AA
        fx = hdulist['FLUX'].data

    # Error array
    if 'ERROR' in hdulist:
        sig = hdulist['ERROR'].data
    else:
        sig = None

    if 'CONTINUUM' in hdulist:
        co = hdulist['CONTINUUM'].data
    else:
        co = None

    xspec1d = XSpectrum1D.from_tuple((wave, fx, sig, co))

    if 'METADATA' in hdulist[0].header:
        # Prepare for JSON (bug fix of sorts)
        metas = hdulist[0].header['METADATA']
        ipos = metas.rfind('}')
        xspec1d.meta.update(json.loads(metas[:ipos+1]))

    return xspec1d


def parse_hdf5(inp, close=True, **kwargs):
    """ Read a spectrum from HDF5 written in XSpectrum1D format
    Expects:  meta, data, units

    Parameters
    ----------
    inp : str or hdf5

    Returns
    -------

    """
    import json
    import h5py
    # Path
    path = kwargs.pop('path', '/')
    # Open
    if isinstance(inp, basestring):
        hdf5 = h5py.File(inp, 'r')
    else:
        hdf5 = inp
    # Data
    data = hdf5[path+'data'].value
    # Meta
    if 'meta' in hdf5[path].keys():
        meta = json.loads(hdf5[path+'meta'].value)
        # Headers
        for jj,heads in enumerate(meta['headers']):
            try:
                meta['headers'][jj] = fits.Header.fromstring(meta['headers'][jj])
            except TypeError:  # dict
                if not isinstance(meta['headers'][jj], dict):
                    raise IOError("Bad meta type")
    else:
        meta = None
    # Units
    units = json.loads(hdf5[path+'units'].value)
    for key,item in units.items():
        if item == 'dimensionless_unit':
            units[key] = u.dimensionless_unscaled
        else:
            units[key] = getattr(u, item)
    # Other arrays
    try:
        sig = data['sig']
    except (NameError, IndexError):
        sig = None
    try:
        co = data['co']
    except (NameError, IndexError):
        co = None
    # Finish
    if close:
        hdf5.close()
    return XSpectrum1D(data['wave'], data['flux'], sig=sig, co=co,
                          meta=meta, units=units, **kwargs)


def parse_DESI_brick(hdulist, select=0):
    """ Read a spectrum from a DESI brick format HDU list

    Parameters
    ----------
    hdulist : FITS HDU list
    select : int, optional
      Spectrum selected. Default is 0

    Returns
    -------
    xspec1d : XSpectrum1D
      Parsed spectrum
    """
    fx = hdulist[0].data
    # Sig
    if hdulist[1].name in ['ERROR', 'SIG']:
        sig = hdulist[1].data
    else:
        ivar = hdulist[1].data
        sig = np.zeros_like(ivar)
        gdi = ivar > 0.
        sig[gdi] = np.sqrt(1./ivar[gdi])
    # Wave
    wave = hdulist[2].data
    wave = give_wv_units(wave)
    if wave.shape != fx.shape:
        wave = np.tile(wave, (fx.shape[0],1))
    # Finish
    xspec1d = XSpectrum1D(wave, fx, sig, select=select)
    return xspec1d


def parse_two_file_format(specfil, hdulist, efil=None):
    """ Parse old two file format (one for flux, another for error).

    Parameters
    ----------
    specfil : str
      Flux filename
    hdulist : FITS HDU list
    efil : str, optional
      Error filename. By default this is inferred from the flux
      filename.

    Returns
    -------
    xspec1d : XSpectrum1D
      Parsed spectrum

    """
    head0 = hdulist[0].header
    # Error
    if efil is None:
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

    # Error file
    if efil is not None:
        efil = os.path.expanduser(efil)
        sig = fits.getdata(efil)
    else:
        sig = None

    #Log-Linear?
    try:
        dc_flag = head0['DC-FLAG']
    except KeyError:
        # The following is necessary for Becker's XShooter output
        cdelt1, dc_flag = get_cdelt_dcflag(head0)

    # Read
    if dc_flag in [0,1]:
        # Data
        fx = hdulist[0].data
        # Generate wave
        wave = setwave(head0)
    else:
        raise ValueError('DC-FLAG has unusual value {:d}'.format(dc_flag))

    # Finish
    xspec1d = XSpectrum1D.from_tuple((wave, fx, sig, None))

    return xspec1d

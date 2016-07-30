"""Module containing the LineLimits object
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import json
import warnings


from astropy import units as u
from astropy.units import Quantity
from astropy import constants as const

class LineLimits(object):
    """ An object for handling the 'limits' of a line

    Parameters
    ----------
    data : `~numpy.ndarray`
        Structured array containing all of the data
        This can be a set of 1D spectra

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object.  e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    Attributes
    ----------
    wvmin : Quantity
      min wavelength (with sig>0.) in selected spectrum
      This may differ from np.min(self.wavelength) if you have not masked the edges
    wvmax : Quantity
      max wavelength (with sig>0.) in selected spectrum
      This may differ from np.min(self.wavelength) if you have not masked the edges
    """

    @classmethod
    def from_file(cls, ifile, **kwargs):
        """ From file

        Parameters
        ----------
        ifile : str
          Filename
        """
        # put import here to avoid circular import with io.py
        #from .io import readspec
        from linetools.spectra import io as tio
        slf = tio.readspec(ifile, **kwargs)
        return slf

    def __init__(self, wave, flux, sig=None, co=None, units=None, select=0,
                 meta=None, verbose=False, masking='edges'):
        """
        Parameters
        ----------
        wave : ndarray  [if 2D, the first axis is nspec]
        flux : ndarray
        sig : ndarray, optional
        units : dict, optional
          Dict containing the units of wavelength, flux
          Required keys are 'wave' and 'flux'
        meta : dict, optional
          Meta data.
          meta['headers'] is a list of input headers (or None's)
        select : int, optional
          Selected Spectrum
        masking: str, optional
          Approach to masking the data using the 'sig' array
          'none'
          'edges' -- Masks all data values with sig <=0 on the 'edge' of each spectrum
             e.g.   sig = [0.,0.,0.,0.2,0.,0.2,0.2,0.,0.] would have the first 3 and last 2 masked
          'all' -- Masks all data values with sig <=0
        """
        # Error checking
        if not isinstance(wave, np.ndarray):
            raise IOError("Input `wave` vector must be an numpy.ndarray.")
        if not isinstance(flux, np.ndarray):
            raise IOError("Input `flux` vector must be an numpy.ndarray.")
        if wave.shape[0] != flux.shape[0]:
            raise IOError("Shape of `flux` and `wave` vectors must be identical.")
        if masking not in ['none', 'edges', 'all']:
            raise IOError("Invalid masking type.")
        #if (masking != 'None') and (sig is None):
        #    warnings.warn("Must input sig array to use masking")

        # Handle many spectra
        if len(wave.shape) == 1:
            self.nspec = 1
            self.totpix = wave.shape[0]
        else:
            self.nspec = wave.shape[0]
            self.totpix = wave.shape[1]
        self.select = select

        if verbose:
            print("We have {:d} spectra with {:d} pixels each.".format(
                self.nspec, self.totpix))

        # Data arrays are always MaskedArray
        self.data = np.ma.empty((self.nspec,), #self.npix),
                               dtype=[(str('wave'), 'float64', (self.totpix)),
                                      (str('flux'), 'float32', (self.totpix)),
                                      (str('sig'),  'float32', (self.totpix)),
                                      (str('co'),   'float32', (self.totpix)),
                                     ])
        self.data['wave'] = np.reshape(wave, (self.nspec, self.totpix))
        self.data['flux'] = np.reshape(flux, (self.nspec, self.totpix))

        if co is not None:
            if wave.shape[0] != co.shape[0]:
                raise IOError("Shape of `wave` and `co` vectors must be identical.")
            self.data['co'] = np.reshape(co, (self.nspec, self.totpix))
        else:
            self.data['co'] = np.nan

        # Need to set sig last for masking
        if sig is not None:
            if wave.shape[0] != sig.shape[0]:
                raise IOError("Shape of `wave` and `sig` vectors must be identical.")
            self.data['sig'] = np.reshape(sig, (self.nspec, self.totpix))
            if masking != 'none':
                for kk in range(self.nspec):
                    gdsigval = np.where(self.data['sig'][kk].data > 0.)[0]
                    badsigval = self.data['sig'][kk].data <= 0.
                    for key in self.data.dtype.names:
                        if masking == 'edges':
                            try:
                                self.data[key][kk][0:gdsigval[0]].mask = True
                            except IndexError:
                                pdb.set_trace()
                            self.data[key][kk][gdsigval[-1]+1:].mask = True
                        elif masking == 'all':
                            self.data[key][kk].mask = badsigval
        else:
            self.data['sig'] = np.nan


        # Units
        if units is not None:
            if not isinstance(units, dict):
                raise IOError("`units` must be a dictionary.")
            valid_keys = ['wave', 'flux']
            for key,item in units.items():
                if key not in valid_keys:
                    raise IOError("{:s} is a wrong key in the `units` dictionary. Valid keys are: {}".format(key, valid_keys))
                assert isinstance(item, UnitBase)
            self.units = units
        else:
            warnings.warn("No unit given to wavelength, assuming Angstroms.")
            self.units = dict(wave=u.AA, flux=u.dimensionless_unscaled)

        # Apply continuum?
        self.normed = False

        # Meta
        if meta is None:
            self.meta = dict(headers=[None]*self.nspec)
        else:
            self.meta = meta
        if 'airvac' not in self.meta.keys():
            self.meta['airvac'] = 'vac'

        # Filename
        self.filename = 'none'

    @property
    def header(self):
        """ Return the header (may be None)
        """
        return self.meta['headers'][self.select]

    @wavelength.setter
    def wavelength(self, value):
        gdp = ~self.data['wave'][self.select].mask
        self.data['wave'][self.select][gdp] = value
        if hasattr(value, 'unit'):
            self.units['wave'] = value.unit


    def vactoair(self):
        """Convert to air-based wavelengths from vacuum

        Returns:
        ----------
        Resultant wavelength array is in AA no matter the input units
        """

    def __dir__(self):
        """ Does something more sensible than what Spectrum1D provides

        Returns
        -------
        dir : list

        """
        return dir(type(self))

    def __repr__(self):
        txt = '<{:s}'.format(self.__class__.__name__)
        # Normalized?
        try:
            if self.normed is True:
                txt = txt + ' (normalized): '
            else:
                txt = txt + ': '
        except AttributeError:
            txt = txt +': '
        # Name
        try:
            txt = txt + 'file={:s},'.format(self.filename)
        except:
            pass
        # nspec, select
        txt = txt + ' nspec={:d},'.format(self.nspec)
        txt = txt + ' select={:d},'.format(self.select)
        # wrest
        txt = txt + ' wvmin={:g}, wvmax={:g}'.format(
            self.wvmin, self.wvmax)
        txt = txt + '>'
        return (txt)

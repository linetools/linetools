"""Module containing the XSpectrum1D class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import json
import warnings


from astropy import units as u
from astropy.units import Quantity, UnitBase
from astropy import constants as const
from astropy.io import fits
#from astropy.nddata import StdDevUncertainty
from astropy.table import QTable, Column, Table

import linetools.utils as liu

from .plotting import get_flux_plotrange
from .utils import meta_to_disk

from ..analysis.interactive_plot import InteractiveCoFit
from ..analysis.continuum import prepare_knots
from ..analysis.continuum import find_continuum
from ..analysis.interp import AkimaSpline
eps = np.finfo(float).eps


class XSpectrum1D(object):
    """ A Class containing 1D spectra and associated methods

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

    @classmethod
    def from_spec1d(cls, spec1d):
        """ Input Spectrum1D

        Avoids error array for now
        """
        slf = cls.from_tuple((spec1d.dispersion, spec1d.flux))
        return slf

    @classmethod
    def from_list(cls, xspecs, **kwargs):
        """ Generate a single XSpectrum1D instance containing an array of
        spectra from a list of individual XSpectrum1D spectra.
        Each spectrum is padded with extra pixels so that the
        wavelength ranges of all spectra are covered.
        Padded pixels are masked.

        Also note that masked pixels in the original data are ignored!

        Uses meta to store headers

        Parameters
        ----------
        xspecs : list

        """
        nspec = len(xspecs)
        # Find max npix
        max_npix = 0
        unit0 = xspecs[0].units
        for xspec in xspecs:
            max_npix = max(max_npix, xspec.npix)
            # Check units
            for key in unit0.keys():
                assert unit0[key] == xspec.units[key]
        # Generate dummy arrays
        wave = np.zeros([nspec, max_npix], dtype='float64')
        flux = np.zeros([nspec, max_npix], dtype='float32')
        if xspec.sig_is_set:
            sig = np.zeros_like(flux)
        else:
            sig = None
        if xspec.co_is_set:
            co = np.zeros_like(flux)
        else:
            co = None
        # Fill
        meta = dict(headers=[])
        for jj,xspec in enumerate(xspecs):
            wave[jj,0:xspec.npix] = xspec.wavelength.value
            flux[jj,0:xspec.npix] = xspec.flux.value
            if xspec.sig_is_set:
                sig[jj,0:xspec.npix] = xspec.sig.value
            if xspec.co_is_set:
                co[jj,0:xspec.npix] = xspec.co.value
            # Meta
            meta['headers'].append(xspec.header)
        # Finish
        return cls(wave,flux,sig=sig,co=co, units=unit0, masking='edges',
                   meta=meta)

    @classmethod
    def from_tuple(cls, ituple, sort=True, **kwargs):
        """Make an XSpectrum1D from a tuple of numpy arrays or Quantity arrays

        Parameters
        ----------
        ituple : (wave,flux), (wave,flux,sig) or (wave,flux,sig,cont)
          ndarray or Quantity array
          If wave is unitless, Angstroms are assumed
          If flux is unitless, it is made dimensionless
          The units for sig and co are taken from flux.
        sort : bool, optional
          Sort by wavelength?
        """
        # Checks
        if not isinstance(ituple,tuple):
            raise IOError("Input tuple only")
        if len(ituple[0].shape) != 1:
            raise IOError("Not ready for multi-dimension.  Sorting will fail..")
        # Parse wavelength
        try:
            wv_unit = ituple[0].unit
        except AttributeError:
            iwave = ituple[0]
            warnings.warn("Assuming wavelength unit is Angstroms")
            wv_unit = u.AA
        else:
            if wv_unit is None:
                wv_unit = u.AA
            try:
                iwave = ituple[0].value
            except AttributeError:
                iwave = ituple[0]

        # Parse flux
        try:
            fx_unit = ituple[1].unit
        except AttributeError:
            iflux = ituple[1]
            fx_unit = u.dimensionless_unscaled
        else:
            if fx_unit is None:
                fx_unit = u.dimensionless_unscaled
            try:
                iflux = ituple[1].value
            except AttributeError:
                iflux = ituple[1]

        # Sort and append None
        ltuple = list(ituple)
        for ii in range(len(ltuple),4):
            ltuple.append(None)
        if sort:
            srt = np.argsort(iwave)
            iwave = iwave[srt]
            iflux = iflux[srt]
            for ii in range(2, len(ltuple)):
                if ltuple[ii] is not None:
                    ltuple[ii] = ituple[ii][srt]

        # Generate
        spec = cls(iwave, iflux, sig=ltuple[2], co=ltuple[3], units=dict(wave=wv_unit, flux=fx_unit), **kwargs)

        # Return
        return spec

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

    def copy(self, select=None):
        """ Copy the spectrum

        Parameters
        ----------
        select : int, optional
          Allow the user to specify a different spectrum on the copy
        """
        # Select
        if select is None:
            select = self.select
        # Key components
        data = self.data.copy()
        units = self.units.copy()
        meta = self.meta.copy()
        #
        if self.sig_is_set:
            sig = data['sig']
        else:
            sig = None
        #
        if self.co_is_set:
            co = data['co']
        else:
            co = None
        new = XSpectrum1D(data['wave'], data['flux'], sig=sig, co=co,
                          units=units, meta=meta, select=select)
        return new

    @property
    def header(self):
        """ Return the header (may be None)
        """
        return self.meta['headers'][self.select]

    @property
    def wavelength(self):
        """ Return the wavelength array with units
        """
        return self.data['wave'][self.select].compressed() * self.units['wave']

    @wavelength.setter
    def wavelength(self, value):
        gdp = ~self.data['wave'][self.select].mask
        self.data['wave'][self.select][gdp] = value
        if hasattr(value, 'unit'):
            self.units['wave'] = value.unit

    @property
    def flux(self):
        """ Return the flux with units
        """
        #flux =  self.data[self.select]['flux'] * self.units['flux']
        flux =  self.data['flux'][self.select].compressed() * self.units['flux']
        if self.normed and self.co_is_set:
            flux /= self.data['co'][self.select].compressed()
        return flux

    @flux.setter
    def flux(self, value):
        gdp = ~self.data['flux'][self.select].mask
        self.data['flux'][self.select][gdp] = value
        if hasattr(value, 'unit'):
            self.units['flux'] = value.unit

    @property
    def sig_is_set(self):
        """ Returns whether the error array is set
        """
        if np.isnan(self.data['sig'][self.select].compressed()[0]):
            return False
        else:
            return True

    @property
    def sig(self):
        """ Return the 1sigma error array.  Adopts the same units as flux
        """
        if not self.sig_is_set:
            warnings.warn("This spectrum does not contain an input error array")
            return np.nan
        #
        #sig = self.data[self.select]['sig'] * self.units['flux']
        sig = self.data['sig'][self.select].compressed() * self.units['flux']
        if self.normed and self.co_is_set:
            sig /= self.data['co'][self.select].compressed()
        return sig

    @sig.setter
    def sig(self, value):
        """ Assumes units are the same as the flux
        """
        gdp = ~self.data['sig'][self.select].mask
        self.data['sig'][self.select][gdp] = value


    @property
    def co_is_set(self):
        """ Returns whether a continuum is defined
        """
        if np.isnan(self.data['co'][self.select].compressed()[0]):
            return False
        else:
            return True

    @property
    def co(self):
        """ Return the continuum.  Assumes the flux units
        """
        if not self.co_is_set:
            warnings.warn("This spectrum does not contain an input continuum array")
            return np.nan
        return self.data['co'][self.select].compressed() * self.units['flux']

    @co.setter
    def co(self, value):
        """ Assumes units are the same as the flux
        """
        gdp = ~self.data['co'][self.select].mask
        self.data['co'][self.select][gdp] = value

    @property
    def npix(self):
        """ Number of *unmasked* pixels """
        try:
            return self._npix
        except AttributeError:
            self.set_diagnostics()
            return self._npix


    @property
    def wvmin(self):
        """Minimum wavelength """
        try:
            return self._wvmin
        except AttributeError:
            self.set_diagnostics()
            return self._wvmin

    @property
    def wvmax(self):
        """Maximum wavelength """
        try:
            return self._wvmax
        except AttributeError:
            self.set_diagnostics()
            return self._wvmax

    def set_diagnostics(self):
        """Generate simple diagnostics on the spectrum.

        As a default, the method cuts on `good' pixels.  Useful for
        plotting, quick comparisons, etc. It currently generates only
        the minimum and maximum wavelengths. Sets attributes `_wvmin`
        and `_wvmax`.
        """
        # Cut on good pixels
        if self.sig_is_set:
            gdpx = self.sig > 0.
        else:
            gdpx = np.array([True] * self.wavelength.value.size)
            #gdpx = np.array([True] * self.data['flux'].size)
        # Fill in attributes
        self._npix = len(self.data['flux'][self.select].compressed())
        self._wvmin = np.min(self.wavelength[gdpx])
        self._wvmax = np.max(self.wavelength[gdpx])

    #  Add noise
    def add_noise(self, seed=None, s2n=None, rstate=None):
        """Add noise to the spectrum.

        Uses the uncertainty array unless s2n is specified.

        Parameters
        ----------
        seed : int, optional
          Seed for the random number generator
        s2n : float, optional
          S/N per pixel for the output spectrum. If None, use the
          uncertainty array.
        rstate : RandomState
          allows regeneration of the same noise

        Returns
        -------
        newspec : A new XSpectrum1D with noise added
        """
        # Seed
        if rstate is None:
            rstate=np.random.RandomState(seed)

        # Random numbers
        #rand = np.random.normal(size=self.npix)
        rand = rstate.randn(self.npix)

        # Modify the flux
        if s2n is not None:
            sig = 1. / s2n
        else:
            sig = self.sig.value
        # Copy
        newspec = self.copy()
        # Deal with mask
        gdp = ~self.data['flux'][self.select].mask
        newspec.data['flux'][self.select][gdp] = self.flux.value + (rand * sig)
        newspec.data['sig'][self.select][gdp] = np.ones_like(self.flux) * sig
        #
        return newspec

    def airtovac(self):
        """ Converts current wavelength array from an assumed air to vacuum wavelength scale

        Returns
        -------
        Resultant wavelength array is in AA no matter the input units
        """
        if self.meta['airvac'] == 'vac':
            warnings.warn("Already in vacuum.  Not applying any correction")
            return
        # Convert to AA
        wavelength = self.wavelength.to(u.AA).value

        # Standard conversion format
        sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
        factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
        factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

        # Convert
        wavelength = wavelength*factor
        # Units
        new_wave = wavelength*u.AA
        # Finish
        self.wavelength = new_wave
        self.meta['airvac'] = 'vac'

    def vactoair(self):
        """Convert to air-based wavelengths from vacuum

        Returns:
        ----------
        Resultant wavelength array is in AA no matter the input units
        """
        # Check
        if self.meta['airvac'] == 'air':
            warnings.warn("Already in air.  Not applying any correction")
            return
        # Convert to AA
        wavelength = self.wavelength.to(u.AA).value

        # Standard conversion format
        sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
        factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
        factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

        # Convert
        wavelength = wavelength/factor
        # Units
        new_wave = wavelength*u.AA
        # Finish
        self.wavelength = new_wave
        self.meta['airvac'] = 'air'

    def constant_sig(self, sigv=0.):
        """Set the uncertainty array to a constant value.

        Parameters
        ----------
        sigv : float, optional
          Scalar sigma value to use.
        """
        self.data['sig'] = sigv

    #  Normalize
    def normalize(self, co, verbose=False, no_check=False):
        """ Normalize the spectrum with an input continuum

        Parameters
        ----------
        co: numpy array
          Continuum.
        verbose: bool [False]
        no_check: bool [False]
          Check size of array?
        """
        if len(co) != self.npix:
            if no_check:
                print('WARNING: Continuum length differs from flux')
                if len(co) > len(self.flux):
                    self.data['co'] = co[0:self.npix]
                else:
                    raise ValueError('normalize: Continuum needs to be longer!')
            else:
                raise ValueError('normalize: Continuum needs to be same length as flux array')
        else:
            self.co = co
        self.normed = True

    def unnormalize(self):
        """ Un-normalize the spectrum to the original flux values

        """
        self.normed = False

    def normalized_spec(self):
        """ Generate and pass back a normalized version of the spectrum
        """
        if self.normed is False:
            raise IOError("Spectrum must be normalized!")
        #
        # Generate
        newspec = XSpectrum1D.from_tuple((self.wavelength, self.flux, self.sig),
                                         meta=self.meta.copy())
        # Return
        return newspec

    def pix_minmax(self, *args):
        """ Find pixel indices given a wavelength or velocity range

        Parameters
        ----------
        Option 1: wvmnx
          * wvmnx: Tuple of 2 floats
            wvmin, wvmax in spectral units

        Option 2: zabs, wrest, vmnx  [not as a tuple or list!]
          * zabs : Absorption redshift
          * wrest : Rest wavelength  (with Units!)
          * vmnx : Tuple/array/list of 2 Quantities
            vmin, vmax in km/s

        Returns
        -------
        gdpix, wvmnx, pixmnx
          * gdpix: Pixel indices satisfying the cut
          * wvmnx: Tuple of the min and max wavelengths
          * pixmnx: Tuple of the min and max pixel indices
        """

        if len(args) == 1:  # Option 1
            wvmnx = args[0]
        elif len(args) == 3:  # Option 2
            from astropy import constants as const
            # args = zabs, wrest, vmnx
            wvmnx = (args[0] + 1) * (args[1] +
                                     (args[1] * args[2] / const.c.to('km/s')))
            wvmnx.to(u.AA)

        # Locate the values
        pixmin = np.argmin(np.fabs(self.wavelength - wvmnx[0]))
        pixmax = np.argmin(np.fabs(self.wavelength - wvmnx[1]))

        gdpix = np.arange(pixmin, pixmax + 1, dtype=int)

        # Fill + Return
        self.sub_pix = gdpix
        return gdpix, wvmnx, (pixmin, pixmax)

    #### ###############################

    # Splice spectrum + Normalize
    # Quick plot
    def plot(self, **kwargs):
        """ Plot the spectrum

        Parameters
        ----------
        show : bool
          If True (the default), then run the matplotlib.pyplot show
          command to display the plot. Disable this if you are running
          from a script and wish to delay showing the plot.
        xlim : tuple of two floats
          The initial x plotting limits (xmin, xmax)
        inline : bool
          Recommended to use if displaying inline in a Notebook
        plot_two : XSpectrum1D
          Plot another spectrum
        scale_two : float
          Scale the 2nd spectrum
        xspec : bool
          Launch XSpecGUI instead

        Other keyword arguments are passed to the matplotlib plot
        command.
        """
        # Launch XSpectrum1D??
        if 'xspec' in kwargs:
            import sys
            from PyQt4 import QtGui
            from linetools.guis.xspecgui import XSpecGui
            app = QtGui.QApplication(sys.argv)
            gui = XSpecGui(self)
            gui.show()
            app.exec_()
            return

        import matplotlib.pyplot as plt
        from ..analysis.interactive_plot import PlotWrapNav
        plt.rcParams['axes.formatter.useoffset'] = False  # avoid scientific notation in axes tick labels

        # Keywords
        nocolor = (False if 'color' in kwargs else True)
        xlim = kwargs.pop('xlim', None)
        inline = kwargs.pop('inline', False)
        xspec2 = kwargs.pop('plot_two', None)
        scale_two = kwargs.pop('scale_two', 1.)

        if inline:
            fig = plt.figure(figsize=(12,8))
        else:
            fig = plt.gcf()
        ax = plt.gca()

        artists = {}
        ax.axhline(0, color='k', lw=0.5)

        show = kwargs.pop('show', True)

        if nocolor:
            kwargs.update(color='0.5')
        artists['fl'] = ax.plot(self.wavelength, self.flux,
                                drawstyle='steps-mid', label='1', **kwargs)[0]

        # Error
        if nocolor:
            kwargs.update(color='g')
        if self.sig_is_set:
            ax.plot(self.wavelength, self.sig, **kwargs)

        # Continuum
        if self.co_is_set and (not self.normed):
            if nocolor:
                kwargs.update(color='r')
            ax.plot(self.wavelength, self.co, **kwargs)

        # Second spectrum
        if xspec2 is not None:
            ax.plot(xspec2.wavelength, xspec2.flux*scale_two, color='blue',
                    label='2')
            legend = ax.legend(loc='upper left', borderpad=0.3,
                            handletextpad=0.3, fontsize='large')

        ax.set_ylim(*get_flux_plotrange(self.flux))

        if xlim is not None:
            xmin, xmax = xlim
        else:
            xmin, xmax = self.wavelength.value[0], self.wavelength.value[-1]

        ax.set_xlim(xmin, xmax)

        # Labels
        ax.set_xlabel('Wavelength ({:s})'.format(self.units['wave'].name), size=16)
        ax.set_ylabel('Flux', size=16)

        if plt.get_backend() == 'MacOSX':
            warnings.warn("""\
            Looks like you're using the MacOSX matplotlib backend. Switch to the TkAgg
            or QtAgg backends to enable all interactive plotting commands.
            """)
        else:
            # Enable xspecplot-style navigation (i/o for zooming, etc).
            # Need to save this as an attribute so it doesn't get
            # garbage-collected.
            self._plotter = PlotWrapNav(
                fig, ax, self.wavelength, self.flux, artists, printhelp=False,
                xlim=(xmin, xmax))

            if show:
                plt.show()

    #  Rebin
    def rebin(self, new_wv, do_sig=False):
        """ Rebin to a new wavelength array

        Uses simple linear interpolation.  The default (and only)
        option conserves counts (and flambda).

        WARNING: Do not trust either edge pixel of the new array.
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

        Returns
        -------
        XSpectrum1D of the rebinned spectrum
        """
        from scipy.interpolate import interp1d
        # Save flux info to avoid unit issues
        funit = self.flux.unit
        flux = self.flux.value

        # Deal with nan
        badf = np.isnan(flux)
        if np.sum(badf) > 0:
            warnings.warn("Ignoring NAN in flux")
        gdf = ~badf
        flux = flux[gdf]

        # Endpoints of original pixels
        npix = len(self.wavelength)
        wvh = (self.wavelength + np.roll(self.wavelength, -1)) / 2.
        wvh[npix - 1] = self.wavelength[npix - 1] + \
                        (self.wavelength[npix - 1] - self.wavelength[npix - 2]) / 2.
        dwv = wvh - np.roll(wvh, 1)
        dwv[0] = 2 * (wvh[0] - self.wavelength[0])
        med_dwv = np.median(dwv.value)

        wvh = wvh[gdf]
        dwv = dwv[gdf]

        # Error
        if do_sig:
            var = self.sig.value**2
            var = var[gdf]
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
        if (bwv[-1] > wvh[-1]):
            newcum[-1] = cumsum[-1]
            newvar[-1] = cumvar[-1]

        # Rebinned flux, var
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
            # Deal with bad pixels
            bad = np.where(var <= 0.)[0]
            for ibad in bad:
                bad_new = np.where(np.abs(new_wv-self.wavelength[ibad]) <
                                   (new_dwv[1:]+dwv[ibad])/2)[0]
                new_sig[bad_new] = 0.
        else:
            new_sig = None

        # update continuum
        if self.co_is_set:
            x, y = self._get_contpoints()
            new_co = self._interp_continuum(x, y, new_wv)
        else:
            new_co = None

        newspec = XSpectrum1D.from_tuple((new_wv, new_fx*funit,
                                          new_sig, new_co),
                                         meta=self.meta.copy())
        # Return
        return newspec

    # Velo array
    def relative_vel(self, wv_obs):
        """ Return a velocity array relative to an input wavelength.

        Should consider adding a velocity array to this Class,
        i.e. self.velo

        Parameters
        ----------
        wv_obs : Quantity
          Wavelength to set the zero of the velocity array.
          Often (1+z)*wrest

        Returns
        -------
        velo : Quantity array (km/s)
        """
        if not isinstance(wv_obs, Quantity):
            raise ValueError('Input wavelength needs to be a Quantity')
        return ((self.wavelength - wv_obs) * const.c / wv_obs).to('km/s')

    #  Box car smooth
    def box_smooth(self, nbox, preserve=False, **kwargs):
        """ Box car smooth the spectrum

        Parameters
        ----------
        nbox: int
          Number of pixels to smooth over
        preserve: bool (False)
          If True, perform a convolution to ensure the new spectrum
          has the same number of pixels as the original.
        **kwargs: dict
          If preserve=True, these keywords are passed on to
          astropy.convoution.convolve

        Returns
        -------
        A new XSpectrum1D instance of the smoothed spectrum
        """
        if preserve:
            from astropy.convolution import convolve, Box1DKernel
            new_fx = convolve(self.flux, Box1DKernel(nbox), **kwargs)
            if self.sig_is_set:
                new_sig = convolve(self.sig, Box1DKernel(nbox), **kwargs)
            else:
                new_sig = None
            new_wv = self.wavelength
        else:
            # Truncate arrays as need be
            npix = len(self.flux)
            try:
                new_npix = npix // nbox  # New division
            except ZeroDivisionError:
                raise ZeroDivisionError('Dividing by zero..')
            orig_pix = np.arange(new_npix * nbox)

            # Rebin (mean)
            new_wv = liu.scipy_rebin(self.wavelength[orig_pix], new_npix)
            new_fx = liu.scipy_rebin(self.flux[orig_pix], new_npix)
            if self.sig_is_set:
                new_sig = liu.scipy_rebin(
                    self.sig[orig_pix], new_npix) / np.sqrt(nbox)
            else:
                new_sig = None

        # Return
        return XSpectrum1D.from_tuple(
            (new_wv, new_fx, new_sig), meta=self.meta.copy())

    # Splice two spectra together
    def gauss_smooth(self, fwhm, **kwargs):
        """ Smooth a spectrum with a Gaussian

        Note that the uncertainty array is not smoothed.

        Parameters
        ----------
        fwhm : float
          FWHM of the Gaussian in pixels (unitless)

        Returns
        -------
        A new XSpectrum1D instance of the smoothed spectrum
        """
        # Import
        from linetools.spectra import convolve as lsc

        # Apply to flux
        new_fx = lsc.convolve_psf(
            self.flux.value, fwhm, **kwargs) * self.flux.unit

        # Get the right sigma
        if self.sig_is_set:
            new_sig = self.sig.value
        else:
            new_sig = None

        # Return
        return XSpectrum1D.from_tuple(
            (self.wavelength, new_fx, new_sig), meta=self.meta.copy())

    def stitch(self, idx=None, scale=1.):
        """ Combine two or more spectra within the .data array
        Simple logic is used to order them by wavelength if the
          order is not specified

        Parameters
        ----------
        idx : list or ndarray
          indices of spectra to stitch and the order to do so
          if None, all of the spectra in the .data array will be combined
            with simple logic using the wavelengths
        scale : float, optional
          Scale factor for flux and error array.

        Returns
        -------
        spec : XSpectrum1D
          The stitched spectrum.
        """
        from linetools.spectra import utils as ltsu
        if idx is None:
            wvmx = []
            for ii in range(self.nspec):
                wvmx.append(np.max(self.data['wave'][ii]))
            # Sort
            idx = np.argsort(np.array(wvmx))
        # Splice the first two
        spec = ltsu.splice_two(self.copy(select=idx[0]),
                               self.copy(select=idx[1]))
        # Loop using the rest
        for kk in range(2,len(idx)):
            spec = ltsu.splice_two(spec.copy(), self.copy(select=idx[kk]))
        # Return
        return spec

    def write(self, outfil, FITS_TABLE=False, **kwargs):
        """  Wrapper for writing
        Parses the extension to choose the file format

        Parameters
        ----------
        outfil : str
          Allowed extensions are
          .fit, .fits -- FITS file; set FITS_TABLE=True to format as a binary FITS Table
          .hdf5 -- HDF5 file
          .ascii -- ASCII
        kwargs

        Returns
        -------

        """
        ext = outfil[outfil.rfind('.')+1:]
        if ext in ['fit','fits']:
            if FITS_TABLE:
                self.write_to_binary_fits_table(outfil, **kwargs)
            else:
                self.write_to_fits(outfil, **kwargs)
        elif ext in ['hdf5']:
            self.write_to_hdf5(outfil, **kwargs)
        elif ext in ['ascii']:
            self.write_to_ascii(outfil, **kwargs)
        else:
            raise IOError("Bad file extension: {:s}".format(ext))

    def write_to_ascii(self, outfil, format='ascii.ecsv', **kwargs):
        """ Write to a text file.

        Parameters
        ----------
        outfil: str
          Filename.
        """
        # Convert to astropy Table
        table = QTable([self.wavelength, self.flux],
                       names=('WAVE', 'FLUX'))
        if self.sig_is_set:
            sigclm = Column(self.sig, name='ERROR')
            table.add_column(sigclm)
        if self.co_is_set:
            coclm = Column(self.co, name='CO')
            table.add_column(coclm)

        # Write
        table.write(outfil, format=format)

    def write_to_fits(self, outfil, select=False, clobber=True, fill_val=0.):
        """ Write to a multi-extension FITS file.

        Writes 2D images for multiple spectra data arrays,
        unless select=True.
        Otherwise writes 1D arrays

        Parameters
        ----------
        outfil : str
          Name of the FITS file
        select : int, optional
          Write only the select spectrum.
          This will always trigger if there is only 1 spectrum
          in the data array
        clobber : bool (True)
          Clobber existing file?
        add_wave : bool (False)
          Force writing of wavelengths as array, instead of using FITS
          header keywords to specify a wcs.
        fill_val : float, optional
          Fill value for masked pixels
        """
        if self.nspec == 1:
            select = True

        # Flux
        if select:
            prihdu = fits.PrimaryHDU(self.data['flux'][self.select].filled(fill_val))
            #prihdu = fits.PrimaryHDU(self.data[self.select]['flux'])
        else:
            prihdu = fits.PrimaryHDU(self.data['flux'].filled(fill_val))
        hdu = fits.HDUList([prihdu])
        prihdu.name = 'FLUX'

        # Wavelength
        if select:
            wvhdu = fits.ImageHDU(self.data['wave'][self.select].filled(fill_val))
        else:
            wvhdu = fits.ImageHDU(self.data['wave'].filled(fill_val))
        wvhdu.name = 'WAVELENGTH'
        hdu.append(wvhdu)

        if self.sig_is_set:
            if select:
                sighdu = fits.ImageHDU(self.data['sig'][self.select].filled(fill_val))
            else:
                sighdu = fits.ImageHDU(self.data['sig'].filled(fill_val))
            sighdu.name = 'ERROR'
            hdu.append(sighdu)

        if self.co_is_set:
            if select:
                cohdu = fits.ImageHDU(self.data['co'][self.select].filled(fill_val))
            else:
                cohdu = fits.ImageHDU(self.data['co'].filled(fill_val))
            cohdu.name = 'CONTINUUM'
            hdu.append(cohdu)

        # Use the header of the selected spectrum
        if self.header is not None:
            hdukeys = list(prihdu.header.keys())
            # Append ones to avoid
            hdukeys = hdukeys + ['BUNIT', 'COMMENT', '', 'NAXIS1', 'NAXIS2', 'HISTORY']
            for key in self.header.keys():
                # Use new ones
                if key in hdukeys:
                    continue
                # Update unused ones
                try:
                    prihdu.header[key] = self.header[key]
                except ValueError:
                    raise ValueError('l.spectra.utils: Bad header key card')
            # History
            if 'HISTORY' in self.header.keys():
                # Strip \n
                tmp = str(self.header['HISTORY']).replace('\n', ' ')
                try:
                    prihdu.header.add_history(str(tmp))
                except ValueError:
                    import pdb
                    pdb.set_trace()

        #
        if self.meta is not None and len(self.meta) > 0:
            prihdu.header['METADATA'] = meta_to_disk(self.meta)

        # Units, etc.
        prihdu.header['NSPEC'] = self.nspec
        prihdu.header['NPIX'] = self.npix
        units = self.units.copy()
        d = liu.jsonify(units)
        # import pdb; pdb.set_trace()
        prihdu.header['UNITS'] = json.dumps(d)

        hdu.writeto(outfil, clobber=clobber)
        print('Wrote spectrum to {:s}'.format(outfil))

    def add_to_hdf5(self, hdf5, path='/', fill_val=0.):
        """ Write the full data array to an already open hdf5 file

        Parameters
        ----------
        hdf5 : h5py.File
          If input, outfil is ignored
        path : str, optional
          Path to the location for writing (useful for using Groups)
        fill_val : float, optional
          Fill value for masked pixels
        """
        # Meta
        if self.meta is not None and len(self.meta) > 0:
            hdf5[path]['meta'] = meta_to_disk(self.meta)
        # Units
        units = self.units.copy()
        d = liu.jsonify(units)
        hdf5[path]['units'] = json.dumps(d)
        # Data with compression
        hdf5.create_dataset(path+'data', data=self.data.filled(fill_val),
                                       chunks=True, compression='gzip')

    def write_to_hdf5(self, outfil, hdf5=None, clobber=True, fill_val=0.):
        """ Write the full data array to an hdf5 file.

        Parameters
        ----------
        outfil : str
          Name of the hdf5
        clobber : bool (True)
          Clobber existing file?
        fill_val : float, optional
          Fill value for masked pixels
        hdf5 : h5py.File, optional
          If input, outfil is ignored
        """
        # Check for h5py
        try:
            import h5py
        except ImportError:
            raise ImportError("You must install h5py to use this method")
        import os
        # Check for file
        if clobber is False:
            if os.path.exists(outfil):
                raise IOError("File exists.  Will only over-write if you set clobber=True")
        # Begin the file
        if hdf5 is None:
            hdf5 = h5py.File(outfil, 'w')
            close = True
        else:
            close = False
        # Meta
        if self.meta is not None and len(self.meta) > 0:
            hdf5['meta'] = meta_to_disk(self.meta)
        # Units
        units = self.units.copy()
        d = liu.jsonify(units)
        hdf5['units'] = json.dumps(d)
        # Data with compression
        hdf5.create_dataset('data', data=self.data.filled(fill_val),
                                       chunks=True, compression='gzip')
        # Finish
        if close:
            hdf5.close()
            print('Wrote spectrum to {:s}'.format(outfil))


    def write_to_binary_fits_table(self, outfil, clobber=True):
        """ Write to a binary FITS table.

        Parameters
        ----------
        outfil : str
          Name of the FITS file
        clobber : bool (True)
          Clobber existing file?
        """
        # TODO
        #  1. Add unit support for wavelength arrays

        # Dummy hdu
        prihdu = fits.PrimaryHDU()
        prihdu.header['NSPEC'] = self.nspec
        prihdu.header['NPIX'] = self.npix

        # TBL HDU
        gdc = [str('wave'),str('flux')]
        if self.sig_is_set:
            gdc += [str('sig')]
        if self.co_is_set:
            gdc += [str('co')]
        tblhdu = fits.BinTableHDU(self.data[gdc].copy(), name='DATA')
        #
        hdulist = fits.HDUList([prihdu, tblhdu])

        # Deal with header
        if hasattr(self, 'head'):
            hdukeys = list(prihdu.header.keys())
            # Append ones to avoid
            hdukeys = hdukeys + ['BUNIT', 'COMMENT', '', 'NAXIS1',
                                 'NAXIS2', 'HISTORY', 'NAXIS', 'END']
            for key in self.header.keys():
                # Use new ones
                if key in hdukeys:
                    continue
                # Update unused ones
                try:
                    prihdu.header[key] = self.header[key]
                except ValueError:
                    raise ValueError('l.spectra.utils: Bad header key card')
            # History
            if 'HISTORY' in self.head.keys():
                # Strip \n
                tmp = str(self.head['HISTORY']).replace('\n', ' ')
                try:
                    prihdu.header.add_history(str(tmp))
                except ValueError:
                    import pdb
                    pdb.set_trace()

        # META
        if self.meta is not None and len(self.meta) > 0:
            prihdu.header['METADATA'] = meta_to_disk(self.meta)

        # Units
        units = self.units.copy()
        d = liu.jsonify(units)
        prihdu.header['UNITS'] = json.dumps(d)

        # Write
        hdulist.writeto(outfil, clobber=clobber)
        print('Wrote spectrum to {:s}'.format(outfil))

    def fit_continuum(self, knots=None, edges=None, wlim=None, dw=10.,
                      kind=None, **kwargs):
        """ Interactively fit a continuum.

        This sets the following attributes
          * spec.co: the new continuum
          * spec.meta['contpoints']: knots defining the continuum

        Use linetools.analysis.interp.AkimaSpline to regenerate the
        continuum from the knots.

        Parameters
        ----------

        knots: list of (x, y) pairs, optional
          A list of spline knots to use for the continuum.
        edges: list of floats, optional
          A list of edges defining wavelength chunks. Spline knots
          will be placed at the centre of these chunks. Default
          is None, which means it will place equally spaced chunks of
          width `dw` (see below).
        wlim : (float, float), optional
          Start and end wavelengths for fitting the continuum. Default is
          None, which fits the entire spectrum.
        dw : float, optional
          The approximate distance between spline knots in
          Angstroms. Default is 10.
        kind : {'QSO', None}, optional
          If not None, generate spline knots using
          linetools.analysis.continuum.find_continuum.
        **kwargs : dict
          Other keyword arguments are passed to
          ~linetools.analysis.continuum.find_continuum.  For
          kind='QSO', allowed keywords are `redshift`, `divmult`,
          `forest_divmult`.

        """
        import matplotlib.pyplot as plt
        if plt.get_backend() == 'MacOSX':
            warnings.warn("""\
            Looks like you're using the MacOSX matplotlib backend. Switch to the TkAgg
            or QtAgg backends to enable all interactive plotting commands.
            """)
            return

        wa = self.wavelength.value
        flux = self.flux.value
        sig = self.sig.value

        anchor = False
        if wlim is None:
            wmin, wmax = wa[0], wa[-1]
        else:
            wmin, wmax = wlim
            if wmax < wmin:
                wmin, wmax = wmax, wmin
            anchor = True

        if kind is not None:
            _, knots = find_continuum(self, kind=kind, **kwargs)
            # only keep knots between wmin and wmax
            knots = [[x, y] for (x, y) in knots if wmin <= x <= wmax]
        else:
            if edges is None:
                nchunks = max(3, (wmax - wmin) / float(dw))
                edges = np.linspace(wmin, wmax, nchunks + 1)

        if knots is None:
            knots, indices, masked = prepare_knots(
                wa, flux, sig, edges)
        else:
            knots = [list(k) for k in knots]

        if not len(knots) > 0:
            raise RuntimeError('Problem generating continuum spline knots.')

        # set the initial continuum for the fitter
        if self.co_is_set:
            co_init = self.co
        else:
            co_init = None
        if co_init is not None:
            x = [k[0] for k in knots]
            ynew = np.interp(x, wa, co_init)
            for i in range(len(knots)):
                knots[i][1] = ynew[i]

        contpoints = [k[:2] for k in knots]
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(11, 7))
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95)
        wrapper = InteractiveCoFit(wa, flux, sig,
                                   contpoints, co=co_init, fig=fig, anchor=anchor)

        # wait until the interactive fitting has finished
        while not wrapper.finished:
            plt.waitforbuttonpress()

        print('Updating continuum.')
        self.co = wrapper.continuum  # this should work with the assignment, even is self.co_is_set is False
        if 'contpoints' not in self.meta:
            self.meta['contpoints'] = []
        self.meta['contpoints'].extend(
            [tuple(pts) for pts in wrapper.contpoints])
        self.meta['contpoints'].sort()

    def _interp_continuum(self, x, y, wv=None):
        """ Interpolate the continuum from spline knots.

        Returns an interpolation of the continuum using knots at
        (x,y) positions. If at least 5 points exists, it uses Akima
        Spline interpolation.  it uses linear interpolation.

        Parameters
        ----------
        x : array, shape(N,)
            The x positions of the knots to interpolate through
            Assumed to be in Angstroms.
        y : array, shape(N,)
            The y positions of the knots to interpolate through
        wv : array, shape(M,); optional
            The wavelength array for final interpolation. If None,
            it uses self.wavelength.value

        Returns
        -------
        co : array, shape(M,)
            Values of the interpolated continuum for each `wv` point

        """
        if wv is None:
            wv = self.wavelength.value

        if len(y) >= 5:
            # need 5 points to define an Akima Spline
            spline = AkimaSpline(x, y)
            co = spline(wv)
        else:
            co = np.interp(wv, x, y)

        return co

    def _get_contpoints(self):
        """Gets the x and y values from the metadata meta['contpoints'] as
        arrays

        Returns
        -------
        x, y : tuple of np.arrays()
            The positions of the contpoint knots in the `x` and `y` axes respectively
        """

        try:
            contpoints = self.meta['contpoints']
        except:
            raise RuntimeError('Meta data `contpoints` does not exist; cannot get contpoints!')
        xy = np.array(contpoints)
        xy = xy.transpose()
        x, y = xy[0], xy[1]
        return x, y

    def perturb_continuum(self, rel_var=0.05, seed=None):
        """ Perturb an existing continuum.

        Adds a perturbation to the continuum by adding noise to the
        contpoints.  This may be useful for estimating the uncertainty
        in the continuum.

        Parameters
        ----------
        rel_var : float, optional
            Relative variation in flux w/r to its original value. This will be
            the standard deviation of a Normal distribution added to
            the contpoints' `y` values. Default is rel_var=0.05 (i.e. 5%).
        seed : int, optional
            The seed of the (pseudo)random number generator.

        Returns
        -------
        Updates self.co with a perturbation.
        Note: To reset to the original continuum use reset_continuum()

        """

        x, y = self._get_contpoints()

        # Get indices in the full spectrum
        # inds = np.array([np.where(np.fabs(x[i]-self.dispersion.value)<0.01)[0][0] for i in range(len(x))])
        # y_err = y_err = self.uncertainty.array[inds]

        # Seed
        if seed is not None:
            np.random.seed(seed=seed)
        #add Normal noise to y points
        y += np.random.normal(0, y * rel_var, len(y))

        #update continuum
        co = self._interp_continuum(x, y, self.wavelength.value)
        self.normalize(co=co)
        #self.co = co

    def reset_continuum(self):
        """Resets the continuum to its original value.

        Returns
        -------
        Updates self.co to its original value.

        Notes
        -----
        This is mainly for use with XSpectrum1D.perturb_continuum
        """

        x, y = self._get_contpoints()

        #update continuum
        co = self._interp_continuum(x, y, self.wavelength.value)
        self.normalize(co=co)
        #self.co = co

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

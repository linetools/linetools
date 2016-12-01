"""Module containing the LineLimits object
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings

from astropy import units as u
from astropy.units import Quantity, UnitConversionError
from astropy import constants as const

from linetools import utils as ltu

ckms = const.c.to('km/s')

class LineLimits(object):
    """ An object for handling the 'limits' of a line

    Properties
    ----------
    z : float
      Redshift
    zlim : tuple of floats
      Redshift limits for a line
      Defined as wave/wrest - 1.
    wvlim : Quantity array
      wavelength limits for the line in observer frame
    vlim : Quantity array
      velocity limits for the line
    """

    @classmethod
    def from_specline(cls, aline, z, zlim):
        """ From AbsLine or Emline

        Parameters
        ----------
        aline : AbsLine
        """
        from ..spectralline import AbsLine, EmLine
        if not isinstance(aline, (AbsLine, EmLine)):
            raise IOError("Input aline must be AbsLine or EmLine")
        #
        slf = cls(aline.wrest, z, zlim)
        return slf

    @classmethod
    def from_dict(cls, idict):
        """ Generate a LineLimits class from an input dict
        keys -- wrest, z, zlim -- are required

        Parameters
        ----------
        idict

        Returns
        -------

        """
        wrest = ltu.convert_quantity_in_dict(idict['wrest'])
        slf = cls(wrest, idict['z'], idict['zlim'])
        # Return
        return slf

    def __init__(self, wrest, z, zlim, **kwargs):
        """
        Parameters
        ----------
        wrest : Quantity
          Rest wavelength of the line.
        z : float
          Redshift of the line.
        zlim : tuple or list
          Redshift limits for a line
          Defined as wave/wrest - 1.
          Ok to have zlim[1]==zlim[0], but then self.is_set() == False
        """
        # Error checking
        if not isinstance(z, (float,int)):
            raise IOError("Input z must be a float")
        if not isinstance(zlim, (tuple,list)):
            raise IOError("Input zlim must be a tuple or list")
        if not isinstance(wrest, Quantity):
            raise IOError("Input wrest must be a quantity")

        # Data
        self._data = {}
        # Set
        self._z = z
        self._wrest = wrest
        self.set(zlim, **kwargs)

    @property
    def z(self):
        """ Return z
        """
        return self._z

    @property
    def wrest(self):
        """ Return wrest
        """
        return self._wrest

    @property
    def zlim(self):
        """ Return zlim
        """
        return self._zlim

    @property
    def wvlim(self):
        """ Return wvlim
        """
        return self._wvlim

    @property
    def vlim(self):
        """ Return vlim
        """
        return self._vlim

    def reset(self):
        """ Update all the values
        """
        #self._data['zlim'] = self._zlim
        self._wvlim = self._wrest*(1+np.array(self._zlim))
        self._vlim = ltu.dv_from_z(self._zlim, self._z)

    def is_set(self):
        """ Query if the limits are set to sensible values
        (i.e. zlim[1]-zlim[0]>0)

        Returns
        -------
        bool

        """
        if self._zlim[1] > self._zlim[0]:
            return True
        else:
            return False

    def set(self, inp, chk_z=False):#, itype='zlim'):
        """ Over-ride = to re-init values

        Parameters
        ----------
        inp : tuple, list, or Quantity array
          * If floats -> zlim : Redshift limits
          * If Quantity array with length units  -> wvlim : Wavelength limits
          * If Quantity array with speed units  -> vlim : Velocity limits
        chk_z : bool, optional
          Demand that zlim bound z
          Often not desired as z can be somewhat arbitrary

        Returns
        -------

        """
        # Checks
        if not isinstance(inp, (tuple, list, Quantity)):
            raise IOError("Input must be tuple, list or Quantity.")
        if np.isclose(self._z, 0.):
            warnings.warn("Redshift=0.  If this is unexpected, set _z and reset limits")

        if isinstance(inp[0], float):  # assume zlim
            self._zlim = inp
        elif isinstance(inp[0], Quantity):  # may be wvlim or vlim
            try:  # assume wvlim
                self._zlim = (inp/self._wrest).decompose().to(
                        u.dimensionless_unscaled).value - 1.
            except UnitConversionError:
                try:  # assume vlim
                    self._zlim = ltu.z_from_dv(inp, self._z)
                except ValueError:
                    raise IOError("Quantity must be length or speed.")
        else:
            raise IOError("Input must be floats or Quantities.")
        # Check
        if chk_z:
            if (self._zlim[0] > self._z) or (self._zlim[1] < self._z):
                #import pdb; pdb.set_trace()
                raise IOError("Invalid input. `zlim` does not bound `z`.")
        # Reset
        self.reset()

    def to_dict(self):
        """ Generate a simple dict of the data

        Returns
        -------
        ldict : dict

        """
        ldict = {}
        ldict['z'] = self.z
        ldict['zlim'] = self.zlim
        for key in ['wvlim', 'vlim', 'wrest']:
            obj = getattr(self, key)
            ldict[key] = dict(value=obj.value,
                              unit=obj.unit.to_string())
        # Return
        return ldict

    def __repr__(self):
        txt = '<{:s}'.format(self.__class__.__name__)
        # wrest
        txt = txt + ' wrest={:g}'.format(self.wrest)
        txt = txt + ' z={:g}'.format(self.z)
        txt = txt + ' zlim={}'.format(self.zlim)
        txt = txt + ' wvlim={}'.format(self.wvlim)
        txt = txt + ' vlim={}'.format(self.vlim)
        txt = txt + '>'
        return (txt)

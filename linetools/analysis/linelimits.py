"""Module containing the LineLimits object
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

from astropy import units as u
from astropy.units import Quantity
from astropy import constants as const

from ..spectralline import AbsLine

ckms = const.c.to('km/s')

class LineLimits(object):
    """ An object for handling the 'limits' of a line

    Properties
    ----------
    zlim : tuple of floats
      Redshift limits for a line
      Defined as wave/wrest - 1.
    wvlim : Quantity array
      wavelength limits for the line
    vlim : Quantity array
      velocity limits for the line in km/s
    """
    @classmethod
    def from_absline(cls, aline, zlim):
        """ From AbsLine

        Parameters
        ----------
        aline : AbsLine
        """
        if not isinstance(aline, float):
            raise IOError("Input aline must be AbsLine")
        #
        slf = cls(aline.wrest, aline.attrib['z'], zlim)
        return slf

    def __init__(self, wrest, z, zlim):
        """
        Parameters
        ----------
        wrest : Quantity
          Rest wavelength of the line.  Should match line.wrest
        z : float
          Redshift of the line.  Should match line.attrib['z']
        zlim : tuple or list
          Redshift limits for a line
          Defined as wave/wrest - 1.
        """
        # Error checking
        if not isinstance(z, float):
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
        self.set(zlim)

    @property
    def zlim(self):
        """ Return zlim
        """
        return self._data['zlim']

    @property
    def wvlim(self):
        """ Return wvlim
        """
        return self._data['wvlim']

    @property
    def vlim(self):
        """ Return vlim
        """
        return self._data['vlim']

    def reset(self):
        """ Update all the values
        """
        self._data['zlim'] = self.zlim
        self._data['wvlim'] = self._wrest*(1+np.array(self.zlim))
        self._data['vlim'] = ckms*((self._data['wvlim']-self._wrest*(1+self._z))/(
            self._wrest*(1+self._z))).decompose()

    def set(self, inp, itype='zlim'):
        """ Over-ride = to re-init values

        Parameters
        ----------
        inp : tuple, list, or Quantity array
          * zlim : Redshift limits
        itype : str, optional
          Input type


        Returns
        -------

        """
        # Checks
        if not isinstance(inp, (tuple, list, Quantity)):
            raise IOError("Input must be tuple, list or Quantity")
        if itype == 'zlim':
            self._data['zlim'] = inp
        else:
            raise IOError("Input type must be zlim, vlim, or wvlim")
        # Reset
        self.reset()

    def __repr__(self):
        txt = '<{:s}'.format(self.__class__.__name__)
        # wrest
        txt = txt + ' wrest={:g}'.format(self._wrest)
        txt = txt + '>'
        return (txt)

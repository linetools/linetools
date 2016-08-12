""" Simple Solar abundance calculations.
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import numbers
import imp

from astropy import constants as const
from astropy.io import ascii
from astropy.utils.misc import isiterable

#from xastropy.xutils import xdebug as xdb
l_path = imp.find_module('linetools')[1]

#
class SolarAbund(object):
    """Class to handle simple Solar Abundance calculations

    Parameters
    ----------
    ref: str, optional
       'Asplund2009' :: Asplund et al. 2009, ARA&A, 47, 481 meteoritic
       table (several photospheric though)
    """
    # Init
    def __init__(self, ref='Asplund2009', verbose=False):

        # Error catching
        if not isinstance(ref, basestring):
            raise TypeError('SolarAbund__init__: Wrong ref type for '
                            'SolarAbund input')
        self.ref = ref

        # Load Data
        print('Loading abundances from {:s}'.format(self.ref))
        self.load_data()
        print('Abundances are relative by number on a '
              'logarithmic scale with H=12') 

    def load_data(self):
        """Grab the Solar Abundance data (in linetools/abund)
        """
        # Data file
        if self.ref == 'Asplund2009':
            dat_file = l_path + '/data/abund/solar_Asplund2009.dat'
            # Read table
            names = ('Elm', 'Abund', 'Z')
            table = ascii.read(dat_file, format='no_header', names=names) 
            # 
        else:
            raise ValueError('Unrecognized reference for SolarAbund: {:s}'.format(self.ref))
        # Save
        self._data = table


    def get_ratio(self, rtio):
        """ Return abundance ratio

        Parameters
        ----------
        rtio : str 
          Element ratio (e.g. 'Si/Fe')        
        """
        # Elements
        elm1, elm2 = rtio.split('/')
        # Abundances
        ab1 = self[elm1]
        ab2 = self[elm2]
        # ratio
        return ab1 - ab2

    def __getitem__(self, k):
        """ Return abundance given an element
 
        Parameters
        ----------
        k : int or str or list/tuple
          * int -- Atomic number (6)
          * str -- Element name (e.g. 'C')

        Returns
        -------
        Abund : float
        """
        # Iterate?
        if isiterable(k) and not isinstance(k, basestring): 
            out_abnd = []
            for ik in k:
                out_abnd.append(self[ik])
            out_abnd = np.array(out_abnd)
            return out_abnd

        if isinstance(k, numbers.Integral): # Atomic number
            mt = np.where(self._data['Z'] == k)[0]
            if len(mt) != 1:
                raise ValueError('Atomic Number not in Table: {:d}'.format(k))
        elif isinstance(k, basestring): # Name
            mt = np.where(self._data['Elm'] == k)[0]
            if len(mt) != 1:
                raise ValueError('Element not in Table: {:s}'.format(k))
        else:
            raise IndexError('Not prepared for this type of input', k)

        # Return
        return self._data['Abund'][mt][0]

    # Printing
    def __repr__(self):
        # Generate sets string
        return '<SolarAbund: {:s}>'.format(self.ref)

"""
Module for LineList Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import os, imp
import copy

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import ascii
from astropy.table import QTable, Table, vstack, Column

#from xastropy.xutils import xdebug as xdb
l_path = imp.find_module('linetools')[1]

#
class SolarAbund(object):
    '''Class to handle simple Solar Abundance calculations

    Parameters:
    ----------
    ref: str, optional
       'Asplund2009' :: AA&RA meteoritic table (several photometric)
    '''
    # Init
    def __init__(self, ref='Asplund2009', verbose=False):

        # Error catching
        if not isinstance(ref, basestring):
            raise TypeError('SolarAbund__init__: Wrong ref type for SolarAbund input')
        self.ref = ref

        # Load Data
        print('Loading abundances from {:s}'.format(self.ref))
        self.load_data()
        print('Abundances are relative, by number on a logarithmic scale with H=12') 

    def load_data(self):
        """Grab the Solar Abundance data (in linetools/abund)
        """
        # Data file
        if self.ref == 'Asplund2009':
            dat_file = l_path+'/data/abund/solar_Asplund2009.dat'
            # Read table
            names=('Elm', 'Abund', 'Z')
            table = ascii.read(dat_file, format='no_header', names=names) 
            # 
        else:
            raise ValueError('Unrecognized reference for SolarAbund: {:s}'.format(self.ref))
        # Save
        self._data = table


    #####
    def __getattr__(self, k):
        ''' Passback an array of the data 
        Parameters:
        ------------
        k: Must be a Column name in the data Table
        '''
        try:
            # First try to access __getattr__ in the parent class.
            # This is needed to avoid an infinite loop 
            out = object.__getattr__(k)
        except AttributeError:
            colm = self._data[k]
            return np.array(self._data[k])
        else:
            return out

    def __getitem__(self, k):
        ''' Passback data as a dict (from the table) for the input line

        Parameters:
        ----------
        k: overloaded
          int -- Atomic number (6)
          str -- Element name (e.g. 'C')
          str -- Element ratio (e.g. 'Si/Fe')

        Returns:
        ----------
        Abund or Abund difference for the ratio
        '''
        if isinstance(k,int): # Atomic number
            mt = np.where(self._data['Z'] == k)[0]
            if len(mt) != 1:
                raise ValueError('Atomic Number not in Table: {:d}'.format(k))
        elif isinstance(k, basestring): # Name
            if '/' in k: # Ratio
                elm1,elm2 = k.split('/')
                ab1 = self[elm1]
                ab2 = self[elm2]
                return ab1-ab2
            else:
                mt = np.where(self._data['Elm'] == k)[0]
                if len(mt) != 1:
                    raise ValueError('Element not in Table: {:s}'.format(k))
        else:
            raise IndexError('Not prepared for this type of input',k)

        # Standard call
        return self._data['Abund'][mt][0]


    # Printing
    def __repr__(self):
        # Generate sets string
        return '[SolarAbund: {:s}]'.format(self.ref)


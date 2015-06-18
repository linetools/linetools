"""
Module for parsing Line List data
  Includes the Dict Definition for the Data
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp

from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.tables import QTable

#from xastropy.xutils import xdebug as xdb
lt_path = imp.find_module('linetools')[1]

#
def line_data():
    ''' Defines the dict for spectral line Data

    Group definition:
    -----------------
        0: None
        1: "All" ISM (intended to be all atomic lines ever observed)
        2: Strong ISM
        4: HI Lyman series
        8: H2
       16: CO
       32: EUV
       64: Galaxy Emission
      128: Galaxy Absorption
      256: AGN
      512: ??
     1024: User1 (Reserved)
     2048: User2 (Reserved)
    '''
    ldict = {
        'name': '',           # Name
        'wrest': 0.*u.AA,     # Rest Wavelength (Quantity)
        'f':  0.,             # Oscillator strength
        'gk': 0.,             # Degeneracy of the upper level
        'Ej': 0.*u.eV,        # Energy of lower level (relative to ground state)
        'Ek': 0.*u.eV,        # Energy of upper level (relative to ground state)
        'Ex': 0./u.cm         # Excitation energy (cm^-1)
        'A': 0./u.s,          # Einstein coefficient
        'gamma': 0.,          # Sum of A
        'nj': 0,              # Orbital level of lower state (or vibrational level)
        'nk': 0,              # Orbital level of upper state (or vibrational level)
        'Jj': 0.,             # Tot ang mom (z projection) of lower state (or rotation level)
        'Jk': 0.,             # Tot ang mom (z projection) of upper state (or rotation level)
        'Z': 0,               # Atomic number (for atoms)
        'ion': 0,             # Ionic state
        'group': 0            # Flag for grouping
        }

#
def read_H2():
    ''' Simple def to read H2 data
    Returns:
    --------
    QTable of H2 lines
    '''
    H2_fil = 
    data = QTable.read(gdfil, guess=False, comment=';')

###########################################

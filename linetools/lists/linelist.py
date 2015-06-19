"""
Module for LineList Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import fits
from astropy.table import QTable, Table, vstack

from xastropy.xutils import xdebug as xdb

#
class LineList(object):
    '''Class to over-load Spectrum1D for new functionality not yet in specutils

    Parameters:
    ----------
    llst_keys: str or list  
      Input to grab line list.  Current options are:
       'ISM'     :: "All" ISM lines
       'Strong'  :: Strong ISM lines
       'HI'      :: HI Lyman series
       'H2'      :: H2 (Lyman-Werner)
       'CO'      :: CO UV band-heads [not yet implemented]
       'EUV'     :: Key EUV lines (for CASBAH project)
       'Gal_E'   :: Galaxy emission lines (HII)
       'Gal_A'   :: Galaxy absorption lines (stellar)
       'AGN'     :: Key AGN lines
    '''
    # Init
    def __init__(self, llst_keys):

        # Error catching
        if type(llst_keys) not in [str,list]:
            raise TypeError('LineList__init__: Wrong type for LineList input')

        # Save
        if type(llst_keys) in [str]:
            self.lists = [llst_keys]
        else:
            self.lists = llst_keys

        # Load Data
        self.load_data()

        # Set lines for use
        self.set_lines()

    # 
    def load_data(self):
        ''' Grab the data for the lines of interest
        '''
        # Import
        from linetools.lists import parse as lilp
        reload(lilp)

        # Define datasets
        dataset = {
            'ism': [None], 
            'molecules': [lilp.read_H2]   # H2 
            }

        # Loop on lists
        for llist in self.lists:
            if llist == 'H2':
                sets = ['molecules']
            else:
                raise ValueError('Not ready for this group')

        full_table = None
        all_func = []
        # Loop on data sets
        for iset in sets:
            # Loop on data sources
            for func in dataset[iset]:
                # Query if read already
                if func not in all_func:
                    # Read
                    table = func()
                    if full_table is None:
                        full_table = table
                    else:
                        full_table = vstack([full_table, table])
                    # Save to avoid repeating
                    all_func.append(func)

        # Save as QTable
        self._fulltable = QTable(full_table)

    #####
    def set_lines(self):
        ''' Parse the lines of interest
        '''
        # Loop on lines
        indices = []
        for llist in self.lists:
            if llist == 'H2':
                gdi = np.where(self._fulltable['mol'] == 'H2')[0]
                if len(gdi) == 0:
                    raise IndexError('set_lines: Found no H2 molecules! Read more data')
                indices.append(gdi)

        # Collate (should grab unique ones!)
        all_idx = np.concatenate( [np.array(itt) for itt in indices] )

        # Parse (consider masking instead)
        self._data = self._fulltable[all_idx]

    #####
    def __getattr__(self,k):
        ''' Passback an array or Column of the data 
        k must be a Column name in the data Table
        '''
        return self._data[k]

    #####
    def __getitem__(self,k, tol=1e-3*u.AA):
        ''' Passback a row of data on the input line
        Paramaters:
        ----------
        k: overloaded
          float,Quantity -- Wavelength
          str -- Name

        Returns:
        ----------
        Astropy Row from the data table
        '''
        if type(k) in [float,Quantity]: # Wavelength
            mt = np.where( np.abs(k-self.wrest) < tol)[0]
        elif type(k) in [str]: # Name
            mt = np.where( k == self.name )[0]
        else:
            raise ValueError('Not prepare for this type')

        # Matches
        if len(mt) == 0:
            print('No such line in the list')
            return None
        elif len(mt) == 1:
            return self._data[mt][0]  # Pass back as a Row not a Table
        else:
            raise ValueError('Multiple lines in the list')


    # Printing
    def __repr__(self):
        # Generate sets string
        for kk,llist in enumerate(self.lists):
            if kk == 0:
                sstr = llist
            else:
                sstr = sstr + ',' + llist
        return '[LineList: {:s}]'.format(sstr)

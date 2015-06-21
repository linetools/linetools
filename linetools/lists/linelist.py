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

#from xastropy.xutils import xdebug as xdb

from linetools.lists import parse as lilp

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
    gd_lines: list, optional
      List of wrest for lines to use (drawn from input linelist)
    '''
    # Init
    def __init__(self, llst_keys, gd_lines=None):

        # Error catching
        if type(llst_keys) not in [str,list,unicode]:
            raise TypeError('LineList__init__: Wrong type for LineList input')

        # Save
        if type(llst_keys) in [str,unicode]:
            self.lists = [llst_keys]
        else:
            self.lists = llst_keys

        # Load Data
        self.load_data()

        # Set lines for use
        self.set_lines(gd_lines=gd_lines)

    # 
    def load_data(self, tol=1e-3*u.AA):
        ''' Grab the data for the lines of interest
        '''
        # Import
        reload(lilp)

        # Define datasets: In order of Priority
        dataset = {
            'ism': [lilp.parse_morton03,lilp.read_verner94], # Verner 1994, Morton 2003 
            'molecules': [lilp.read_H2]   # H2 
            }

        # Loop on lists
        sets = []
        for llist in self.lists:
            if str(llist) == 'H2':
                sets.append('molecules')
            elif str(llist) == 'ISM':
                sets.append('ism')
            elif str(llist) == 'Strong':
                sets.append('ism')
            else:
                import pdb
                pdb.set_trace()
                raise ValueError('load_data: Not ready for this: {:s}'.format(llist))

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
                        # Unique values
                        wrest = full_table['wrest']
                        newi = []
                        for jj,row in enumerate(QTable(table)): # QTable for units
                            mt = np.where(np.abs(row['wrest']-wrest) < tol)[0]
                            if len(mt) == 0:
                                newi.append(jj)
                        # Append new ones (can't stack QTables yet)
                        full_table = vstack([full_table, table[newi]])
                    # Save to avoid repeating
                    all_func.append(func)

        # Save as QTable
        self._fulltable = QTable(full_table)

    #####
    def set_lines(self, verbose=True, gd_lines=None):
        ''' Parse the lines of interest
        Parameters:
        -------------
        gd_lines: list, optional
          List of wrest for lines to use (drawn from input linelist)
          Should be unitless, i.e. not Quantity
        '''

        indices = []
        set_flags = []
        if gd_lines is None:  # Default list
            # Loop on lines
            for llist in self.lists:
                if llist == 'H2':
                    gdi = np.where(self._fulltable['mol'] == 'H2')[0]
                    if len(gdi) == 0:
                        raise IndexError(
                            'set_lines: Found no H2 molecules! Read more data')
                    indices.append(gdi)
                elif llist == 'ISM':
                    set_flags.append('fISM')
                elif llist == 'Strong':
                    set_flags.append('fSI')
                else:
                    raise ValueError('set_lines: Not ready for this: {:s}'.format(llist))
        else: # Input subset of lines
            wrest = self._fulltable['wrest'].value # Assuming Anstroms
            for gdlin in gd_lines:
                mt = np.where( 
                    np.abs(gdlin-wrest) < 1e-3 )[0]
                if len(mt) > 0:
                    indices.append(mt)
                else:
                    if verbose:
                        print('set_lines: Did not find {:g} in data Tables'.format(gdlin))

        # Deal with Defined sets
        if len(set_flags) > 0:
            # Read standard file
            set_data = lilp.read_sets()
            # Speed up
            wrest = self._fulltable['wrest'].value # Assuming Anstroms
            for sflag in set_flags:
                gdset = np.where(set_data[sflag] == 1)[0]
                # Match to wavelengths
                for igd in gdset:
                    mt = np.where( 
                        np.abs(set_data[igd]['wrest']-wrest) < 1e-3 )[0]
                    if len(mt) > 0:
                        for imt in mt:
                            # Over-ride name!
                            self._fulltable[imt]['name'] = set_data[igd]['name']
                            #if set_data[igd]['name'] == 'DI 1215':
                            #    xdb.set_trace()
                        indices.append(mt)
                    else:
                        if verbose:
                            print('set_lines: Did not find {:s} in data Tables'.format(
                                set_data[igd]['name']))

        # Collate (should grab unique ones!)
        all_idx = np.concatenate( [np.array(itt) for itt in indices] )

        # Parse (consider masking instead)
        self._data = self._fulltable[all_idx]

    #####
    def __getattr__(self,k):
        ''' Passback an array or Column of the data 
        k must be a Column name in the data Table
        '''
        # Deal with QTable
        colm = self._data[k]
        if type(colm[0]) is Quantity:
            return self._data[k]
        else:
            return np.array(self._data[k])

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
        if type(k) in [float]: # Wavelength, assuming Ang
            mt = np.where( np.abs(k*u.AA-self.wrest) < tol)[0]
        elif type(k) in [Quantity]: # Wavelength
            mt = np.where( np.abs(k-self.wrest) < tol)[0]
        elif type(k) in [str,unicode]: # Name
            mt = np.where( str(k) == self.name )[0]
        else:
            raise ValueError('Not prepare for this type')

        # Matches
        if len(mt) == 0:
            print('No such line in the list')
            return None
        elif len(mt) == 1:
            return self._data[mt][0]  # Pass back as a Row not a Table
        else:
            raise ValueError(
                '{:s}: Multiple lines in the list'.format(self.__class__))


    # Printing
    def __repr__(self):
        # Generate sets string
        for kk,llist in enumerate(self.lists):
            if kk == 0:
                sstr = llist
            else:
                sstr = sstr + ',' + llist
        return '[LineList: {:s}]'.format(sstr)

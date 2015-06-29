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
       'CO'      :: CO UV band-heads
       ---- NOT IMPLEMENTED YET -----
       'EUV'     :: Key EUV lines (for CASBAH project)
       'Gal_E'   :: Galaxy emission lines (HII)
       'Gal_A'   :: Galaxy absorption lines (stellar)
       'AGN'     :: Key AGN lines
    gd_lines: list, optional
      List of wrest for lines to use (drawn from input linelist)
    verbose: bool, optional
      Give info galore if True
    '''
    # Init
    def __init__(self, llst_keys, gd_lines=None, verbose=False):

        # Error catching
        if not isinstance(llst_keys,(str,list,unicode)):
            raise TypeError('LineList__init__: Wrong type for LineList input')

        # Save
        if isinstance(llst_keys, basestring):
            self.lists = [llst_keys]
        else:
            self.lists = llst_keys

        # Take closest line?
        self.closest = False

        # Load Data
        self.load_data()

        # Set lines for use
        self.set_lines(gd_lines=gd_lines, verbose=verbose)

    # 
    def load_data(self, tol=1e-3*u.AA):
        ''' Grab the data for the lines of interest
        '''
        # Import
        reload(lilp)

        # Define datasets: In order of Priority
        dataset = {
            'ism': [lilp.parse_morton03,lilp.parse_morton00, 
                lilp.read_verner94], # Morton 2003, Morton 00, Verner 94 
            'molecules': [lilp.read_H2,lilp.read_CO]   # H2 (Abrigail), CO (JXP)
            }

        # Loop on lists
        sets = []
        flag_fval = False # Update f-values?
        flag_wrest = False # Update wavelengths?
        for llist in self.lists:
            if str(llist) == 'H2':
                sets.append('molecules')
            elif str(llist) == 'CO':
                sets.append('molecules')
            elif str(llist) == 'ISM':
                sets.append('ism')
                flag_fval = True
                flag_wrest = True
            elif str(llist) == 'Strong':
                sets.append('ism')
                flag_fval = True
                flag_wrest = True
            elif str(llist) == 'HI':
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
                            mt = np.abs(row['wrest']-wrest) < tol
                            if mt.sum() == 0:
                                newi.append(jj)
                        # Append new ones (can't stack QTables yet)
                        full_table = vstack([full_table, table[newi]])
                    # Save to avoid repeating
                    all_func.append(func)

        # Save as QTable
        self._fulltable = QTable(full_table)

        # Update wavelength values
        if flag_wrest:
            lilp.update_wrest(self._fulltable)

        # Update f-values (Howk00)
        if flag_fval:
            lilp.update_fval(self._fulltable)
        #import pdb
        #pdb.set_trace()

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

        # Default list
        if gd_lines is None:  
            # Loop on lines
            for llist in self.lists:
                if llist in ['H2','CO']:
                    gdi = np.where(self._fulltable['mol'] == llist)[0]
                    if len(gdi) == 0:
                        raise IndexError(
                            'set_lines: Found no {:s} molecules! Read more data'.format(llist))
                    indices.append(gdi)
                elif llist == 'ISM':
                    set_flags.append('fISM')
                elif llist == 'Strong':
                    set_flags.append('fSI')
                elif llist == 'HI':
                    set_flags.append('fHI')
                else:
                    raise ValueError('set_lines: Not ready for this: {:s}'.format(llist))
        else: # Input subset of lines
            wrest = self._fulltable['wrest'].value # Assuming Anstroms
            for gdlin in gd_lines:
                mt = np.where( 
                    np.abs(gdlin-wrest) < 1e-4 )[0]
                if len(mt) == 1:
                    indices.append(mt)
                elif len(mt) > 1:
                    import pdb
                    pdb.set_trace()
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
                        np.abs(set_data[igd]['wrest']-wrest) < 1e-4 )[0]
                    if len(mt) == 1:
                        for imt in mt:
                            # Over-ride name!
                            self._fulltable[imt]['name'] = set_data[igd]['name']
                            #if set_data[igd]['name'] == 'DI 1215':
                            #    xdb.set_trace()
                        indices.append(mt)
                    elif len(mt) > 1:
                        print('wrest = {:g}'.format(set_data[igd]['wrest']))
                        import pdb
                        pdb.set_trace()
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
        if isinstance(colm[0], Quantity):
            return self._data[k]
        else:
            return np.array(self._data[k])

    #####
    def __getitem__(self,k, tol=1e-3*u.AA):
        ''' Passback data as a dict (from the table) for the input line

        Parameters:
        ----------
        k: overloaded
          float,Quantity -- Wavelength
          str -- Name

        Returns:
        ----------
        Dict (from row in the data table)
        '''
        if isinstance(k,(float,Quantity)): # Wavelength
            if isinstance(k,float): # Assuming Ang
                inwv = k*u.AA
            else:
                inwv = k
            mt = np.where( np.abs(inwv-self.wrest) < tol)[0]
        elif isinstance(k, basestring): # Name
            mt = np.where( str(k) == self.name )[0]
        else:
            raise ValueError('Not prepared for this type')

        # No Match?
        if len(mt) == 0:
            # Take closest??
            if self.closest and (not isinstance(k, basestring)):
                mt = [np.argmin(np.abs(inwv-self.wrest))]
                print('WARNING: Using {:g} for your input {:g}'.format(self.wrest[mt[0]], 
                    inwv))
            else:
                print('No such line in the list')
                return None

        # Now we have something
        if len(mt) == 1:
            return dict(zip(self._data.dtype.names,self._data[mt][0]))  # Pass back as dict
            #return self._data[mt][0]  # Pass back as a Row not a Table
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

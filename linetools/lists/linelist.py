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
       'ISM'     :: "All" ISM lines (can be overwhelming!)
       'Strong'  :: Strong ISM lines
       'HI'      :: HI Lyman series
       'H2'      :: H2 (Lyman-Werner)
       'CO'      :: CO UV band-heads
       'EUV'     :: EUV lines (for CASBAH project)
       ---- NOT IMPLEMENTED YET -----
       'Gal_E'   :: Galaxy emission lines (HII)
       'Gal_A'   :: Galaxy absorption lines (stellar)
       'AGN'     :: Key AGN lines
    gd_lines: list, optional
      List of wrest for lines to use (drawn from input linelist)
      Needs to be Quantity
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
        """Grab the data for the lines of interest
        """
        # Import
        reload(lilp)

        # Define datasets: In order of Priority
        dataset = {
            'ism': [lilp.parse_morton03,lilp.parse_morton00, 
                lilp.read_verner94, lilp.read_euv], # Morton 2003, Morton 00, Verner 94, Verner 96 [soon]
            'molecules': [lilp.read_H2,lilp.read_CO],   # H2 (Abrigail), CO (JXP)
            'euv': [lilp.read_euv]   # EUV lines (by hand for now; soon to be Verner96)
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
            elif str(llist) == 'EUV':
                sets.append('ism')
                sets.append('euv')
                flag_fval = True
                flag_wrest = True
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
        import warnings

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
                elif llist == 'EUV':
                    set_flags.append('fEUV')
                else:
                    raise ValueError('set_lines: Not ready for this: {:s}'.format(llist))
        else: # Input subset of lines
            wrest = self._fulltable['wrest'].value # Assuming Angstroms
            for gdlin in gd_lines:
                mt = np.where( 
                    np.abs(gdlin.value-wrest) < 1e-4 )[0]
                if len(mt) == 1:
                    indices.append(mt)
                elif len(mt) > 1:
                    raise ValueError('Need unique entries!')
                else:
                    if verbose:
                        print('set_lines: Did not find {:g} in data Tables'.format(gdlin))

        # Deal with Defined sets
        #import pdb
        #pdb.set_trace()
        if len(set_flags) > 0:
            # Read standard file
            set_data = lilp.read_sets()
            # Speed up
            wrest = self._fulltable['wrest'].value # Assuming Angstroms
            for sflag in set_flags:
                gdset = np.where(set_data[sflag] == 1)[0]
                # Match to wavelengths
                for igd in gdset:
                    mt = np.where( 
                        np.abs(set_data[igd]['wrest']-wrest) < 9e-5 )[0]
                    if len(mt) == 1:
                        for imt in mt:
                            # Over-ride name!
                            self._fulltable[imt]['name'] = set_data[igd]['name']
                            #if set_data[igd]['name'] == 'DI 1215':
                            #    xdb.set_trace()
                        indices.append(mt[0])
                    elif len(mt) > 1:
                        #
                        wmsg = 'WARNING: Multiple lines with wrest={:g}'.format(
                            set_data[igd]['wrest'])
                        warnings.warn(wmsg)
                        warnings.warn('Taking the first entry. Maybe use higher precision.')
                        indices.append(mt[0])
                    else:
                        if verbose:
                            print('set_lines: Did not find {:s} in data Tables'.format(
                                set_data[igd]['name']))

        # Collate (should grab unique ones!)
        #all_idx = np.unique( np.concatenate( [np.array(itt) for itt in indices] ) )
        all_idx = np.unique( np.array(indices) )

        # Parse and sort (consider masking instead)
        tmp_tab = self._fulltable[all_idx]
        tmp_tab.sort('wrest')
        #
        self._data = tmp_tab

    def unknown_line(self):
        """Returns a dictionary of line properties set to an unknown
        line. Currently using the default value from ."""     
        ldict , _ = lilp.line_data()
        ldict['name'] = 'unknown'
        return ldict

    def all_transitions(self,line):
        """For a given single line transition, this function returns a
        all transitions of the ion containing such single line found in 
        the linelist.

        Parameters:
        ----------
        line: string
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII')
        
            [Note: when string contains spaces it only considers the first
             part of it, so 'HI' and 'HI 1215' and 'HI 1025' are all equivalent]
            [Note: to retrieve an unknown line use string 'unknown']

        Returns:
        ----------
        dict (if only 1 transition found) or Table (if > 1 transitions are found)

        """

        if isinstance(line, basestring): # Name
            line = line.split(' ')[0] # keep only the first part of input name
        else:
            raise ValueError('Not prepared for this type')

        if line == 'unknown':
            return self.unknown_line()
        else:
            Z = None
            data = self._data
            for row in data: #is this loop avoidable?
                name = row['name']
                #keep only the first part of name in linelist too
                name = name.split(' ')[0]
                if name == line:
                    Z = row['Z'] #atomic number
                    ie = row['ion'] #ionization estate
                    break
            if Z is not None:
                tbl = self.__getitem__((Z,ie))
                # For hydrogen/deuterium this contains deuterium/hydrogen; 
                # so let's get rid of them
                if (line == 'HI') or (line == 'DI'):
                    names = np.array(tbl['name'])
                    cond = np.array([l.startswith(line) for l in names])
                    tbl = tbl[cond]
                if len(tbl) > 1:
                    return tbl
                else: #this whould be always len(tbl)==1 because Z is not None
                    name = tbl['name'][0]
                    return self.__getitem__(name)
            else:
                raise ValueError('Line {} not found in the linelist'.format(line))

            
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
          float,Quantity -- Wavelength (e.g. 1215.6700)
          str -- Name (e.g. 'CII 1334')
          tuple -- Zion, e.g. (6,2)
          [Note: to retrieve an unknown line use string 'unknown']

        Returns:
        ----------
        dict (from row in the data table) or Table (tuple)
        '''
        if isinstance(k,(float,Quantity)): # Wavelength
            if isinstance(k,float): # Assuming Ang
                inwv = k*u.AA
            else:
                inwv = k
            mt = np.where( np.abs(inwv-self.wrest) < tol)[0]
        elif isinstance(k, basestring): # Name
            if k == 'unknown':
                return self.unknown_line()
            else:
                mt = np.where(str(k) == self.name)[0]
        elif isinstance(k, tuple): # Zion
            mt = (self._data['Z'] == k[0]) & (self._data['ion'] == k[1])
        else:
            raise ValueError('Not prepared for this type',k)

        # No Match?
        if len(mt) == 0:
            # Take closest??
            if self.closest and (not isinstance(k, basestring)):
                mt = [np.argmin(np.abs(inwv-self.wrest))]
                print('WARNING: Using {:.4f} for your input {:.4f}'.format(self.wrest[mt[0]], 
                    inwv))
            else:
                print('No such line in the list', k)
                return None

        # Now we have something
        if len(mt) == 1:
            return dict(zip(self._data.dtype.names,self._data[mt][0]))  # Pass back as dict
            #return self._data[mt][0]  # Pass back as a Row not a Table
        elif isinstance(k,tuple):
            return self._data[mt]
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

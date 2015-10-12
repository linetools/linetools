"""
Module for LineList Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import copy

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import fits
from astropy.table import QTable, Table, vstack, Column

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
    subset: list, optional
      List of subset of lines to use (drawn from input linelist)
      Needs to be Quantity or str (e.g. [1215.6700*u.AA] or ['HI 1215'])
    sort_subset: bool, optional
      Sort the subset? [False]
    verbose: bool, optional
      Give info galore if True
    '''
    # Init
    def __init__(self, llst_keys, subset=None, verbose=False,
        sort_subset=False):

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

        # Set lines for use (from defined LineList)
        self.set_lines(verbose=verbose)

        # Subset of lines for use
        if subset is not None:
            self.subset_lines(subset, verbose=verbose, sort=sort_subset)

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
            'hi': [lilp.parse_morton03],
            'molecules': [lilp.read_H2,lilp.read_CO],   # H2 (Abrigail), CO (JXP)
            'euv': [lilp.read_euv]   # EUV lines (by hand for now; soon to be Verner96)
            }

        # Loop on lists
        sets = []
        flag_fval = False # Update f-values?
        flag_wrest = False # Update wavelengths?
        flag_gamma = True # Update gamma values (recommended)
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
                sets.append('hi')
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

        # Update gamma-values (Mainly HI)
        if flag_gamma:
            lilp.update_gamma(self._fulltable)

    #####
    def set_lines(self, verbose=True):#, gd_lines=None):
        ''' Parse the lines of interest
        '''
        import warnings

        indices = []
        set_flags = []

        # Default list
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

    def subset_lines(self, subset, sort=False, reset_data=False, verbose=False):
        '''
        Select a user-specific subset of the lines from the LineList for usage

        Parameters:
        -------------
        subset: list, optional
          List of wrest for lines to use (drawn from input linelist)
          Quantity or str
        reset_data: bool, optional
          Reset self._data QTable based on the original list at the 
          initialization (i.e. the default list). This is useful for 
          changing subsets of lines without the need to initialize a 
          different LineList() object. [False]
        sort: bool, optional 
          Sort this subset? [False]

        '''
        
        # Reset _data (useful for changing subsets)
        if reset_data:
            self.set_lines(verbose=False)

        indices = []
        if isinstance(subset[0],(float,Quantity)): # wrest
            wrest = self._data['wrest'].value # Assuming Angstroms
            for gdlin in subset:
                mt = np.where( 
                    np.abs(gdlin.value-wrest) < 1e-4 )[0]
                if len(mt) == 1:
                    indices.append(mt)
                elif len(mt) > 1:
                    raise ValueError('Need unique entries!')
                else:
                    if verbose:
                        print('subset_lines: Did not find {:g} in data Tables'.format(gdlin))
        elif isinstance(subset[0],(basestring)): # Names
            names = np.array(self._data['name'])
            for gdlin in subset:
                mt = np.where(str(gdlin)==names)[0]
                if len(mt) == 1:
                    indices.append(mt[0])
                elif len(mt) > 1:
                    raise ValueError('Should have been only one line with name {:s}!'.format(str(gdlin)))
                    #warnings.warn('Found more than one line for {:s}'.format(str(gdlin)))
                    #warnings.warn('Taking the first one from Ref={:s}'.format(
                    #        self._data['Ref'][mt[0]]))
                    #indices.append(mt[0])
                    #raise ValueError('Need unique name entries!')
                else:
                    if verbose:
                        print('subset_lines: Did not find {:s} in data Tables'.format(gdlin))
            #import pdb
            #pdb.set_trace()
        else:
            raise ValueError('Not ready for this type of gd_lines')
        # Sort
        tmp = self._data[np.array(indices)]
        if sort:
            tmp.sort('wrest')
        # Finish
        self._data = tmp

    def unknown_line(self):
        """Returns a dictionary of line properties set to an unknown
        line. Currently using the default value from linetools.lists.parse()."""     
        ldict , _ = lilp.line_data()
        ldict['name'] = 'unknown'
        return ldict

    def all_transitions(self,line):
        """For a given single line transition, this function returns
        all transitions of the ion containing such single line found in 
        the linelist.

        Parameters
        ----------
        line: str or Quantity
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII', 1215.6700*u.AA)
        
            [Note: when string contains spaces it only considers the first
             part of it, so 'HI' and 'HI 1215' and 'HI 1025' are all equivalent]
            [Note: to retrieve an unknown line use string 'unknown']

        Returns
        -------
        dict (if only 1 transition found) or QTable (if > 1 transitions are found)

        """

        if isinstance(line, basestring): # Name
            line = line.split(' ')[0] # keep only the first part of input name
        elif isinstance(line, Quantity): # Rest wavelength (with units)
            data = self.__getitem__(line)
            return self.all_transitions(data['name'])
        else:
            raise ValueError('Not prepared for this type')

        if line == 'unknown':
            return self.unknown_line()
        else:
            Z = None
            for row in self._data: #is this loop avoidable?
                name = row['name']
                #keep only the first part of name in linelist too
                name = name.split(' ')[0]
                if name == line:
                    Z = row['Z'] #atomic number
                    ie = row['ion'] #ionization estate
                    Ej = row['Ej'] #Energy of lower level
                    break
            if Z is not None:
                tbl = self.__getitem__((Z,ie))
                # Make sure the lower energy level is the same too
                cond = tbl['Ej'] == Ej
                tbl = tbl[cond]
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

    def strongest_transitions(self,line, wvlims, n_max=3,verbose=False):
        """For a given single line transition, this function returns
        the n_max strongest transitions of the ion species found in
        the linelist, within the wavelenght range wlims.

        Parameters
        ----------
        line: str or Quantity
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII',
            1215.6700*u.AA) [Note: when string contains spaces it only
            considers the first part of it, so 'HI' and 'HI 1215' and
            'HI 1025' are all equivalent] [Note: to retrieve an
            unknown line use string 'unknown']
        wvlims : tuple of Quantity, or Quantity tuple
            Wavelength range, e.g. wvlims=(1100*u.AA, 3200*u.AA) or
            wvlims=(1100, 3200)*u.AA
        n_max : int or None
            Maximum number of transitions to retrieve; if n_max=None
            it retrieves all of them

        Returns
        -------
        None (if no transitions are found), dict (if only 1 transition
        found), or QTable (if > 1 transitions are found)

        """    

        #Check correct format
        if isinstance(wvlims, tuple): # tuple
            if all(isinstance(wvlim,Quantity) for wvlim in wvlims):
                pass
            else:
                raise SyntaxError('Elements of wvlims have to be of class Quantity; correct format please')
        elif isinstance(wvlims, Quantity): # or quantity
            pass
        else:
            raise SyntaxError('wvlims has to be tuple or Quantity')
        if len(wvlims) != 2:
            raise SyntaxError('wlims has to be of shape (2,); Please correct format')
        if wvlims[0] >= wvlims[1]:
            raise ValueError('Minimum limit `wlims[0]` is not smaller than maximum limit `wlims[1]`; please correct')
        if isinstance(n_max,int):
            if n_max < 1:
                return None
        elif (n_max is not None):
            raise SyntaxError('n_max must be integer or None')

        data = self.all_transitions(line)
        # condition to be within wvrange
        cond = (data['wrest'] >= wvlims[0]) & (data['wrest'] <= wvlims[1])
        if np.sum(cond) == 0:
            if verbose:
                print('[strongest_transitions] Warning: no transitions found within wvlims; returning None')
            return None
        elif isinstance(data,dict): #Only 1 case from a dict format
            return data
        elif np.sum(cond) == 1: #only 1 case from a QTable format
            name = data[cond]['name'][0]
            return self.__getitem__(name)
        else:
            #remove transitions out of range
            data = data[cond]
            #sort by strength defined as wrest * fosc
            strength = data['wrest'] * data['f']
            sorted_inds = np.argsort(strength)
            #reverse sorted indices, so strongest get first
            sorted_inds = sorted_inds[::-1]
            #sort using sorted_inds
            data = data[sorted_inds]            
            #keep only the first n_max or less
            if n_max is not None:
                data = data[:n_max]            
            if len(data) == 1: #Only 1 case from a QTable format; return a dictionary
                name = data['name'][0]
                return self.__getitem__(name)
            else:
                return data

    def available_transitions(self, wvlims, n_max=None,n_max_tuple=None, min_strength=1.):
        """For a given wavelength range, wvlims=(wv_min,wv_max), this
        function retrieves the n_max_tuple strongest transitions per
        each ion species in the LineList available at such a
        wavelength range and having strength larger than min_strength.
        Strength is defined as log10(wrest*fosc*abundance). The output
        is sorted by strength of the strongest available transition
        per ion species.

        Parameters
        ----------
        wvlims : tuple of Quantity
            Wavelength range, e.g. wvlims=(1100*u.AA, 3200*u.AA)
        n_max : int, optional
            Maximum number of transitions retrieved when given,
            otherwise recover all of them
        n_max_tuple : int, optional
            Maximum number of transitions in a given ion species to
            retrieve. e.g., if Lyman series are all available, it will
            retrieve only up to Lyman gamma if
            n_max_tuple=3. Otherwise it returns all of them
        min_strength : float, optional
            Minimum strength calculated from log10(wrest * fosc *
            abundance) In thin space HI 1215 has 14.7.

        Returns
        -------
        dict (if only 1 transition found) or QTable (if > 1
        transitions are found) or None (if no transition is found)
        """
        if all((isinstance(n,int) or (n is None)) for n in [n_max,n_max_tuple]):
            if (n_max is not None) and (n_max < 1):
                return None
        else:
            raise SyntaxError('Both n_max and n_max_tuple must be integers when given!')
        if isinstance(min_strength,float) or isinstance(min_strength,int):
            pass
        else:
            raise SyntaxError('min_strength must be a float value')

        # Identify unique ion_names (e.g. HI, CIV, CIII)
        #unique_ion_names = list(set([name.split(' ')[0] for name in self._data['name']]))
        #unique_ion_names = np.array(unique_ion_names)
        unique_ion_names = np.unique([name.split(' ')[0] for name in self._data['name']])

        #obtain the strongest transition of a given unique ion species
        ion_name = []
        strength = []
        for ion in unique_ion_names: #This loop is necesary to have a non trivial but convinient order in the final output
            aux = self.strongest_transitions(ion,wvlims,n_max=1) #only the strongest
            if aux is not None:
                if isinstance(aux,dict):#this should always be True given n_max=1
                    name = aux['name']
                else:
                    name = aux['name'][0]
                abundance = get_abundance(name)[0]
                ion_name += [name]
                strength += [np.log10(aux['wrest'].value * aux['f']) + abundance]
        if len(ion_name)==0:
            #no matches
            return None

        #create Table
        unique = Table()
        unique.add_column(Column(data=ion_name,name='name'))
        unique.add_column(Column(data=strength,name='strength'))

        #get rid of those below the min_strength threshold
        cond = unique['strength'] >= min_strength
        unique = unique[cond]
        if len(unique) < 1:
            return None

        #sort by strength
        unique.sort(['strength'])
        unique.reverse() #Table unique is now sorted by strength, with only 
                         #1 entry per ion species

        #Create output data adding up to n_max_tuple per ion species
        for i,row in enumerate(unique):
            name = row['name']
            aux = self.strongest_transitions(name, wvlims, n_max=n_max_tuple)
            #need to deal with dict vs QTable format now
            if isinstance(aux,dict):
                aux = self.from_dict_to_qtable(aux)
            if i == 0:
                output = Table(aux) #convert to table because QTable does not like vstack
            else:
                output = vstack([output,Table(aux)]) #vstack is only supported for Table()
        #if len==1 return dict
        if len(output) == 1:
            name = output['name'][0]
            return self.__getitem__(name)
        else: #n_max>1
            if n_max>1:
                output = output[:n_max]
            if len(output) == 1: #return dictionary
                name = output['name'][0]
                return self.__getitem__(name)
            else:
                return QTable(output)

    def from_dict_to_qtable(self,a):
        """Converts dictionary `a` to its QTable version"""
        if isinstance(a,dict):
            pass
        else:
            raise SyntaxError('Input has to be a dictionary')
        
        keys = self._data.keys()

        #Create a QTable with same shape as self._data
        tab = QTable(data=self._data[0])
        #re-write the value elements
        for key in keys:
            tab[0][key] = a[key]
        return tab
    
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
        dict (from row in the data table if only 1 line is found) or QTable (tuple
          when more than 1 lines are found)
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


def get_abundance(transitions_names): #please remove this function later 
    """Temporary function to obtain an array of abundances from 
    a given array of ion names. The abundance scheme should be 
    implemented in a better way!!! This is only a temporary function"""

    #Create a dictionary of abundances [temporary, until abundances 
    # is implemented at a higher level somewhere]

    if isinstance(transitions_names,basestring):
        transitions_names = [transitions_names]

    transitions_names = np.array(transitions_names)
    

    abundance = {
    'H': 12.00,
    'He': 10.9,
    'Li': 1.05,
    'Be': 1.38,
    'B': 2.70,
    'C': 8.43,
    'N': 7.83,
    'O': 8.69,
    'F': 4.56,
    'Ne': 7.93,
    'Na': 6.24,
    'Mg': 7.60,
    'Al': 6.45,
    'Si': 7.51,
    'P': 5.41,
    'S': 7.12,
    'Cl': 5.50,
    'Ar': 6.40,#jump to Fe
    'Fe': 7.50
    } #assume the rest are very small for now; see below

    abund = []
    #create abundance array
    for name in transitions_names: 
        #keep only the element
        ion = name.split(' ')[0]
        #get atom name
        if ion[1].islower():
            atom = ion[:2]
        else:
            atom = ion[0]
        #check whether atom is in abundance dictionary
        if atom in abundance.keys():
            abund += [abundance[atom]]
        else: #if not in key, use a very low abundance [Temporary!!!]
            abund += [1.00]
    return np.array(abund)

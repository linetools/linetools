""" Contains the LineList class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.table import QTable, Table, vstack, Column, MaskedColumn
import warnings

import imp

lt_path = imp.find_module('linetools')[1]

# from xastropy.xutils import xdebug as xdb

CACHE = {'full_table': {}, 'data': {}}

from linetools.lists import parse as lilp

# TODO
# Do something about transitions that are both in Galaxy and ISM
# (e.g. MgII). Currently, priority is given to the first one loaded
# but that may not be the best approach...


class LineList(object):
    """
    This Class is designed to organize and handle information about
    atomic and/or molecular transition lines (e.g. HI Lya, CIV 1548,
    Hydrogen Balmer series, etc.) observed in a variety of
    astrophysical environments. It is currently implemented for
    absorption lines, but we expect to also include common emission
    lines in the near future.

    Parameters
    ----------
    llst_key : str
        Input to grab line list.  Current options are:
        * 'ISM'     :: "All" ISM lines (can be overwhelming!)
        * 'Strong'  :: Strong ISM lines
        * 'HI'      :: HI Lyman series
        * 'H2'      :: H2 (Lyman-Werner)
        * 'CO'      :: CO UV band-heads
        * 'EUV'     :: EUV lines (for CASBAH project)
        * 'Galaxy'  :: Lines typically identified in galaxy spectra
        * 'AGN'     :: Key AGN lines (to be implemented)

    verbose : bool, optional
        Give info galore if True

    use_ISM_table : bool [default True]
        Developer use only. Read from a stored fits table with all ISM
        transitions, rather than from the original source files.

    sort_key : list of str
        Name(s) of the key(s) to sort the LineList by. Default is [`wrest`].
    """

    # Init
    def __init__(self, llst_key, verbose=False, closest=False, set_lines=True,
                 use_ISM_table=True, use_cache=True, sort_by=['wrest']):

        # Error catching
        if not isinstance(llst_key, basestring):
            raise TypeError('LineList__init__: Wrong type for LineList input')

        # Save
        self.list = llst_key

        # Take closest line?
        self.closest = closest
        self.verbose = verbose

        if not use_ISM_table or llst_key not in ('ISM', 'HI', 'Strong', 'EUV'):
            # Load Data
            self.load_data(use_cache=use_cache)

        if set_lines:
            # Set lines for use (from defined LineList)
            # This sets self._data
            self.set_lines(use_ISM_table=use_ISM_table, verbose=verbose,
                           use_cache=use_cache)
        # Memoize
        self.memoize = {}  # To speed up multiple calls

        # set strength (using default values for now)
        self._set_extra_columns_to_datatable()

        # sort the LineList
        self.sort(sort_by)
        self.sort_by = sort_by

    def load_data(self, use_ISM_table=True, tol=1e-3 * u.AA, use_cache=True):
        """Grab the data for the lines of interest
        """
        global CACHE
        key = self.list, tol
        if use_cache and key in CACHE['full_table']:
            self._fulltable = CACHE['full_table'][key]
            return

        # Define datasets: In order of Priority
        dataset = {
            'ism': [lilp.parse_morton03, lilp.parse_morton00, lilp.parse_verner96,
                    lilp.read_verner94, lilp.read_euv],  # Morton 2003, Morton 00, Verner 96, Verner 94
            'hi': [lilp.parse_morton03],
            # H2 (Abrigail), CO (JXP)
            'molecules': [lilp.read_H2, lilp.read_CO],
            'euv': [lilp.read_euv],  # EUV lines (by hand for now; soon to be Verner96)
            'galaxy': [lilp.read_forbidden, lilp.read_recomb, lilp.read_galabs],
        }

        sets = []
        flag_fval = False  # Update f-values?
        flag_wrest = False  # Update wavelengths?
        flag_gamma = True  # Update gamma values (recommended)

        if self.list in ('H2', 'CO'):
            sets.append('molecules')
        elif self.list in ('ISM', 'Strong', 'EUV'):
            flag_fval = True
            flag_wrest = True
            sets.append('ism')
            if self.list == 'EUV':
                sets.append('euv')
        elif self.list == 'HI':
            sets.append('hi')
        elif self.list == 'Galaxy':
            sets.append('galaxy')
        else:
            # import pdb
            # pdb.set_trace()
            raise ValueError(
                'load_data: Not ready for this: {:s}'.format(self.list))

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
                    # make it a QTable
                    table = QTable(table)

                    if full_table is None:
                        full_table = table
                    else:
                        # Unique values
                        wrest = full_table['wrest']
                        newi = []
                        for jj, row in enumerate(table):
                            try:
                                mt = np.abs(row['wrest'] - wrest) < tol
                            except:
                                import pdb
                                pdb.set_trace()
                            if mt.sum() == 0:
                                newi.append(jj)
                        # Append new ones as Tables (can't stack QTables yet)
                        full_table = vstack([Table(full_table), Table(table[newi])])
                    # Save to avoid repeating
                    all_func.append(func)

        # Save as QTable
        self._fulltable = QTable(full_table)

        # Update wavelength values
        if flag_wrest:
            lilp.update_wrest(self._fulltable)

        # Update f-values (Howk00)
        if flag_fval:
            self._fulltable = lilp.update_fval(self._fulltable)

        # Update gamma-values (Mainly HI)
        if flag_gamma:
            lilp.update_gamma(self._fulltable)

        CACHE['full_table'][key] = self._fulltable

    #####
    def set_lines(self, verbose=True, use_cache=True, use_ISM_table=True):
        """ Parse the lines of interest

        Parameters
        ----------
        verbose : bool, optional
        use_cache : bool, optional
          cache the linelist for faster repeat performance
        use_ISM_table : bool, optional
          For speed, use a saved ISM table instead of reading from original source files.
        """

        global CACHE
        key = self.list
        if use_cache and key in CACHE['data']:
            self._data = CACHE['data'][key]
            return
        elif use_ISM_table and self.list in ('ISM', 'Strong', 'EUV', 'HI'):
            data = QTable(Table.read(lt_path + '/data/lines/ISM_table.fits'))
            if self.list != 'ISM':
                cond = data['is_'  + self.list]
                self._data = data[cond]
            else:
                self._data = data
            CACHE['data'][key] = self._data
            return

        indices = []
        set_flags = []

        # Default list
        # Loop on lines
        if self.list in ['H2', 'CO']:
            gdi = np.where(self._fulltable['mol'] == self.list)[0]
            if len(gdi) == 0:
                raise IndexError(
                    'set_lines: Found no {:s} molecules! Read more data'.format(self.list))
            indices.append(gdi)
        elif self.list == 'ISM':
            set_flags.append('fISM')
        elif self.list == 'Strong':
            set_flags.append('fSI')
        elif self.list == 'HI':
            set_flags.append('fHI')
        elif self.list == 'EUV':
            set_flags.append('fEUV')
        elif self.list == 'Galaxy':
            set_flags.append('fgE')
            set_flags.append('fgA')
        else:
            raise ValueError(
                'set_lines: Not ready for this: {:s}'.format(self.list))

        # Deal with Defined sets
        # import pdb
        # pdb.set_trace()
        if len(set_flags) > 0:
            # Read standard file
            set_data = lilp.read_sets()
            # Speed up
            wrest = self._fulltable['wrest'].value  # Assuming Angstroms
            for sflag in set_flags:
                gdset = np.where(set_data[sflag] == 1)[0]
                # Match to wavelengths
                for igd in gdset:
                    mt = np.where(
                        np.abs(set_data[igd]['wrest'] - wrest) < 9e-5)[0]
                    if len(mt) == 1:
                        for imt in mt:
                            # Over-ride name!
                            self._fulltable[imt][
                                'name'] = set_data[igd]['name']
                            # if set_data[igd]['name'] == 'DI 1215':
                            #    xdb.set_trace()
                        indices.append(mt[0])
                    elif len(mt) > 1:
                        #
                        wmsg = 'WARNING: Multiple lines with wrest={:g}'.format(
                            set_data[igd]['wrest'])
                        warnings.warn(wmsg)
                        warnings.warn(
                            'Taking the first entry. Maybe use higher precision.')
                        indices.append(mt[0])
                    else:
                        if verbose:
                            print('set_lines: Did not find {:s} in data Tables'.format(
                                set_data[igd]['name']))

        # Collate (should grab unique ones!)
        all_idx = np.unique(np.array(indices))

        # Parse and sort (consider masking instead)
        tmp_tab = self._fulltable[all_idx]
        tmp_tab.sort('wrest')

        #
        self._data = tmp_tab
        CACHE['data'][key] = self._data

    def _set_extra_columns_to_datatable(self, abundance_type='solar', ion_correction='none'):
        """Sets new convenient columns to the self._data QTable.
        These include:
            - `ion_name` : HI, CIII, CIV, etc
            - `log(w*f)` : np.log10(wrest * fosc)  # in np.log10(AA)
            - `abundance` : given by [`none`, `solar`]
            - `ion_correction` : [`none`]
            - `rel_strength` : log(w*f) + abundance + ion_correction

        This function is only implemented for the following lists: HI, ISM, EUV, Strong

        """
        good_linelists = ['HI', 'ISM', 'EUV', 'Strong']
        if self.list not in good_linelists:
            warnings.warn('Not implemented: will not set relative strength for LineList: {}.'.format(self.list))
            return

        # Set ion_name column
        ion_name = [name.split(' ')[0] for name in self._data['name']]
        self._data['ion_name'] = ion_name

        # Set QM strength as MaskedColumn (self._data['f'] is MaskedColumn)
        qm_strength = self._data['f'] * (self._data['wrest'].to('AA').value)
        qm_strength.name = 'qm_strength'
        self._data['log(w*f)'] = np.log10(qm_strength)
        # mask out potential nans
        cond = np.isnan(self._data['log(w*f)'])
        self._data['log(w*f)'].mask = np.where(cond, True, self._data['log(w*f)'].mask)

        # Set Abundance
        available_abundance_types = ['none', 'solar']
        if abundance_type not in available_abundance_types:
            raise ValueError('_set_extra_columns_to_datatable: `abundance type` '
                             'has to be either: {}'.format(available_abundance_types))
        if abundance_type == 'none':
            self._data['abundance'] = np.ones(len(self._data))
        elif abundance_type == 'solar':
            from linetools.abund.solar import SolarAbund
            solar = SolarAbund()
            # abund will be masked array as default (all masked out) in
            # case an element is not in SolarAbund()
            abund = np.ma.masked_array(np.zeros(len(self._data)), mask=True)
            for ii in range(len(abund)):
                ion_name = self._data['name'][ii]
                ion_Z = self._data['Z'][ii]
                if ion_name.startswith('DI'): # Deuterium
                    abund[ii] = solar['D']
                    abund.mask[ii] = False  # unmask
                else:
                    try:
                        abund[ii] = solar[ion_Z]
                        abund.mask[ii] = False  # unmask
                    except ValueError:
                        pass  # these correspond to elements with no abundance given by solar()
                              # and so they remain masked
            self._data['abundance'] = abund

        # Set ionization correction
        available_ion_corrections = ['none']
        if ion_correction not in available_ion_corrections:
            raise ValueError('_set_extra_columns_to_datatable: `ion_correction` '
                             'has to be either: {}'.format(available_ion_corrections))
        if ion_correction == 'none':
            self._data['ion_correction'] = np.zeros(len(self._data))  # is in log10 scale, so 0 means no-correction

        # Set relative strength in log10 scale
        self._data['rel_strength'] = self._data['log(w*f)'] + self._data['abundance'] + self._data['ion_correction']

    def sort(self, keys, reverse=False):
        """Sort the LineList according to a given key or keys.

        Parameters
        ----------
        key : str or list of str
            The main key(s) to sort the LineList by
            (e.g. 'wrest' or ['Z', 'log(w*f)']).
        reverse : bool, optional
            If True, the sorting is reversed

        Note: this is a wrapper to astropy.table.table.sort()
        """
        # define the sorting key(s) as list
        if isinstance(keys, (str, basestring)):
            keys = [keys]

        # sort
        self._data.sort(keys)

        # reverse?
        if reverse:
            self._data.reverse()

        # update sort_by
        self.sort_by = keys

    def subset_lines(self, subset, reset_data=False, verbose=False, sort_by=['wrest']):
        """ Select a user-specific subset of the lines from the LineList

        Parameters
        ----------
        subset : list or np.ndarray (of Quantity or str)
            List of wrest or names for lines to use (drawn from current
            LineList object) e.g. (['HI 1215', 'CIV 1548'] or
            [1215.67 * u.AA, 1548.195 * u.AA])
        reset_data : bool, optional
            Reset self._data QTable based on the original list at the
            initialization(i.e. the default list). This is useful for
            changing subsets of lines without the need to initialize a
            different LineList() object. Default is False.
        sort_by : list of str, optional
            Key to sort the LineList by

        Returns
        -------
        A new LineList object with the subset of transitions

        """

        # Check the right format
        if not isinstance(subset, (list, np.ndarray)):
            raise ValueError('subset_lines: the input subset must be a list!')

        # Reset _data (useful for changing subsets)
        if reset_data:
            self.set_lines(verbose=False)

        indices = []
        if isinstance(subset[0], (float, Quantity)):  # wrest
            wrest = self._data['wrest'].to('AA').value  # In Angstroms
            for gdlin in subset:
                mt = np.where(
                    np.abs(gdlin.to('AA').value - wrest) < 1e-4)[0]
                if len(mt) == 1:
                    indices.append(mt[0])
                    # import pdb; pdb.set_trace()
                elif len(mt) > 1:
                    raise ValueError('There are multiple matches for line {:g} {:s}!'.format(gdlin.value, gdlin.unit))
                else:
                    if verbose:
                        print(
                            'subset_lines: Did not find {:g} in data Tables'.format(gdlin))
        elif isinstance(subset[0], (basestring)):  # Names
            names = np.array(self._data['name'])
            for gdlin in subset:
                mt = np.where(str(gdlin) == names)[0]
                if len(mt) == 1:
                    indices.append(mt[0])
                    # import pdb; pdb.set_trace()
                elif len(mt) > 1:
                    raise ValueError(
                        'Should have been only one line with name {:s}!'.format(str(gdlin)))
                else:
                    if verbose:
                        print(
                            'subset_lines: Did not find {:s} in data Tables'.format(gdlin))
        else:
            raise ValueError('Not ready for this `subset` type yet.')

        # Return LineList object
        new = LineList(self.list, closest=self.closest, set_lines=False,
                       verbose=self.verbose, sort_by=sort_by)
        new._data = tmp
        new.sort(sort_by)
        return new

    def unknown_line(self):
        """Returns a dictionary of line properties set to an unknown
        line.

        Currently using the default value from
        linetools.lists.parse().
        """
        ldict, _ = lilp.line_data()
        ldict['name'] = 'unknown'
        return ldict

    def all_transitions(self, line):
        """ Get all the transitions corresponding to a line.

        For a given single line transition, this function returns
        all transitions from the LineList containing that line.

        Parameters
        ----------
        line : str or Quantity
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII',
            1215.6700*u.AA). When string contains spaces it only
            considers the first part of it, so 'HI' and 'HI 1215' and
            'HI 1025' are all equivalent. To retrieve an unknown line
            use string 'unknown'.

        Returns
        -------
        dict or QTable
            dict if only 1 transition found, otherwise QTable.

        """
        if isinstance(line, basestring):  # Name
            line = line.split(' ')[0]  # keep only the first part of input name
        elif isinstance(line, Quantity):  # Rest wavelength (with units)
            data = self.__getitem__(line)
            return self.all_transitions(data['name'])
        else:
            raise ValueError('Not prepared for this type')

        if line == 'unknown':
            return self.unknown_line()
        else:
            Z = None
            for row in self._data:  # is this loop avoidable?
                name = row['name']
                # keep only the first part of name in linelist too
                name = name.split(' ')[0]
                if name == line:
                    Z = row['Z']  # atomic number
                    ie = row['ion']  # ionization estate
                    Ej = row['Ej']  # Energy of lower level
                    break
            if Z is not None:
                tbl = self.__getitem__((Z, ie))
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
                else:  # this whould be always len(tbl)==1 because Z is not None
                    name = tbl['name'][0]
                    return self.__getitem__(name)
            else:
                raise ValueError(
                    'Line {} not found in the linelist'.format(line))

    def strongest_transitions(self, line, wvlims, n_max=3, verbose=False):
        """ Find the strongest transition for an ion

        For a given single line transition, this function returns
        the n_max strongest transitions of the ion species found in
        the linelist, within the wavelength range wlims.

        Parameters
        ----------
        line : str or Quantity
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII',
            1215.6700 * u.AA)[Note: when string contains spaces it only
            considers the first part of it, so 'HI' and 'HI 1215' and
            'HI 1025' are all equivalent][Note: to retrieve an
            unknown line use string 'unknown']
        wvlims : tuple of Quantity, or Quantity tuple
            Wavelength range, e.g. wvlims = (1100 * u.AA, 3200 * u.AA) or
            wvlims = (1100, 3200) * u.AA
        n_max : int or None, optional
            Maximum number of transitions to retrieve; if n_max = None
            it retrieves all of them

        Returns
        -------
        None (if no transitions are found), dict (if only 1 transition
        found), or QTable (if > 1 transitions are found)

        """

        # Check correct format
        if isinstance(wvlims, tuple):  # tuple
            if all(isinstance(wvlim, Quantity) for wvlim in wvlims):
                pass
            else:
                raise SyntaxError(
                    'Elements of wvlims have to be of class Quantity; correct format please')
        elif isinstance(wvlims, Quantity):  # or quantity
            pass
        else:
            raise SyntaxError('wvlims has to be tuple or Quantity')
        if len(wvlims) != 2:
            raise SyntaxError(
                'wlims has to be of shape (2,); Please correct format')
        if wvlims[0] >= wvlims[1]:
            raise ValueError(
                'Minimum limit `wlims[0]` is not smaller than maximum limit `wlims[1]`; please correct')
        if isinstance(n_max, int):
            if n_max < 1:
                return None
        elif (n_max is not None):
            raise SyntaxError('n_max must be integer or None')

        data = self.all_transitions(line)
        # condition to be within wvrange
        cond = (data['wrest'] >= wvlims[0]) & (data['wrest'] <= wvlims[1])
        if np.sum(cond) == 0:
            if verbose:
                print(
                    '[strongest_transitions] Warning: no transitions found within wvlims; returning None')
            return None
        elif isinstance(data, dict):  # Only 1 case from a dict format
            return data
        elif np.sum(cond) == 1:  # only 1 case from a QTable format
            name = data[cond]['name'][0]
            return self.__getitem__(name)
        else:
            # remove transitions out of range
            data = data[cond]
            # sort by strength defined as wrest * fosc
            strength = data['wrest'] * data['f']
            sorted_inds = np.argsort(strength)
            # reverse sorted indices, so strongest get first
            sorted_inds = sorted_inds[::-1]
            # sort using sorted_inds
            data = data[sorted_inds]
            # keep only the first n_max or less
            if n_max is not None:
                data = data[:n_max]
            if len(data) == 1:  # Only 1 case from a QTable format; return a dictionary
                name = data['name'][0]
                return self.__getitem__(name)
            else:
                return data

    def available_transitions(self, wvlims, n_max=None, n_max_tuple=None, min_strength=1.):
        """ Find the strongest transitions in a wavelength interval.

        For a given wavelength range, wvlims = (wv_min, wv_max), this
        function retrieves the n_max_tuple strongest transitions per
        each ion species in the LineList available at such a
        wavelength range and having strength larger than min_strength.
        Strength is defined as log10(wrest * fosc * abundance). The output
        is sorted by strength of the strongest available transition
        per ion species.

        Parameters
        ----------
        wvlims : tuple of Quantity
            Wavelength range, e.g. wvlims = (1100 * u.AA, 3200 * u.AA)
        n_max : int, optional
            Maximum number of transitions retrieved when given,
            otherwise recover all of them
        n_max_tuple : int, optional
            Maximum number of transitions in a given ion species to
            retrieve. e.g., if Lyman series are all available, it will
            retrieve only up to Lyman gamma if
            n_max_tuple = 3. Otherwise it returns all of them
        min_strength : float, optional
            Minimum strength calculated from log10(wrest * fosc *
            abundance) In this way HI 1215 has 14.7 by definition.

        Returns
        -------
        dict (if only 1 transition found) or QTable (if > 1
        transitions are found) or None (if no transition is found)
        """
        # Init
        from linetools.abund.solar import SolarAbund
        from linetools.abund import ions as laions
        solar = SolarAbund()

        if all((isinstance(n, int) or (n is None)) for n in [n_max, n_max_tuple]):
            if (n_max is not None) and (n_max < 1):
                return None
        else:
            raise SyntaxError(
                'Both n_max and n_max_tuple must be integers when given!')
        if isinstance(min_strength, (float,int)):
            pass
        else:
            raise SyntaxError('min_strength must be a float value')

        # Identify unique ion_names (e.g. HI, CIV, CIII)
        # unique_ion_names = list(set([name.split(' ')[0] for name in self._data['name']]))
        # unique_ion_names = np.array(unique_ion_names)
        unique_ion_names = np.unique(
            [name.split(' ')[0] for name in self._data['name']])

        # obtain the strongest transition of a given unique ion species
        ion_name = []
        strength = []
        for ion in unique_ion_names:  # This loop is necessary to have a non trivial but convinient order in the final output
            # Abundance
            Zion = laions.name_ion(ion)
            if ion == 'DI':
                abundance = 12. - 4.8  # Approximate for Deuterium
            else:
                abundance = solar[Zion[0]]

            aux = self.strongest_transitions(
                ion, wvlims, n_max=1)  # only the strongest
            if aux is not None:
                if isinstance(aux, dict):  # this should always be True given n_max=1
                    name = aux['name']
                else:
                    name = aux['name'][0]
                ion_name += [name]
                strength += [np.log10(aux['wrest'].value *
                                      aux['f']) + abundance]
        if len(ion_name) == 0:
            # no matches
            return None

        # create Table
        unique = Table()
        unique.add_column(Column(data=ion_name, name='name'))
        unique.add_column(Column(data=strength, name='strength'))

        # get rid of those below the min_strength threshold
        cond = unique['strength'] >= min_strength
        unique = unique[cond]
        if len(unique) < 1:
            return None

        # sort by strength
        unique.sort(['strength'])
        unique.reverse()  # Table unique is now sorted by strength, with only
        # 1 entry per ion species

        # Create output data adding up to n_max_tuple per ion species
        for i, row in enumerate(unique):
            name = row['name']
            aux = self.strongest_transitions(name, wvlims, n_max=n_max_tuple)
            # need to deal with dict vs QTable format now
            if isinstance(aux, dict):
                aux = self.from_dict_to_qtable(aux)
            if i == 0:
                # convert to Table because QTable does not like vstack
                output = Table(aux)
            else:
                # vstack is only supported for Table()
                output = vstack([output, Table(aux)])
        # if len==1 return dict
        if len(output) == 1:
            name = output['name'][0]
            return self.__getitem__(name)
        else:  # n_max>1
            if (n_max is not None) and (n_max > 1):
                output = output[:n_max]
            if len(output) == 1:  # return dictionary
                name = output['name'][0]
                return self.__getitem__(name)
            else:
                return QTable(output)

    def from_dict_to_qtable(self, a):
        """Converts dictionary `a` to its QTable version.
        """

        if isinstance(a, dict):
            pass
        else:
            raise SyntaxError('Input has to be a dictionary')

        keys = self._data.keys()

        # Create a QTable with same shape as self._data
        tab = QTable(data=self._data[0])
        # re-write the value elements
        for key in keys:
            tab[0][key] = a[key]
        return tab

    #####
    def __getattr__(self, k):
        """ Passback an array or Column of the data

        k must be a Column name in the data Table
        """
        # Deal with QTable
        try:
            # First try to access __getitem__ in the parent class.
            # This is needed to avoid an infinite loop which happens
            # when trying to assign self._fulldata to the cache
            # dictionary.
            out = object.__getattr__(k)
        except AttributeError:
            colm = self._data[k]
            if isinstance(colm[0], Quantity):
                return self._data[k]
            else:
                return np.array(self._data[k])
        else:
            return out

    def __getitem__(self, k, tol=1e-3*u.AA):
        """ Passback data as a dict (from the table) for the input line

        Parameters
        ----------
        k: overloaded
          * float, Quantity -- Wavelength (e.g. 1215.6700)
          * str -- Name (e.g. 'CII 1334')
          * tuple -- Zion, e.g. (6,2)
          [Note: to retrieve an unknown line use string 'unknown']

        Returns
        -------
        dict (from row in the data table if only 1 line is found) or
          QTable (tuple when more than 1 lines are found)
        """
        try:
            tmp = self.memoize[k].copy()
        except KeyError:
            if isinstance(k, (float, Quantity)):  # Wavelength
                if isinstance(k, float):  # Assuming Ang
                    inwv = k * u.AA
                else:
                    inwv = k
                mt = np.where(np.abs(inwv - self.wrest) < tol)[0]
            elif isinstance(k, basestring):  # Name
                if k == 'unknown':
                    return self.unknown_line()
                else:
                    mt = np.where(self.name == str(k))[0]
            elif isinstance(k, tuple):  # Zion
                mt = (self._data['Z'] == k[0]) & (self._data['ion'] == k[1])
            else:
                raise ValueError('Not prepared for this type', k)

            # No Match?
            if len(mt) == 0:
                # Take closest??
                if self.closest and (not isinstance(k, basestring)):
                    mt = [np.argmin(np.abs(inwv - self.wrest))]
                    if self.verbose:
                        print('WARNING: Using {:.4f} for your input {:.4f}'.format(self.wrest[mt[0]],
                                                                               inwv))
                else:
                    if self.verbose:
                        print('No such line in the list', k)
                    return None

            # Now we have something
            if len(mt) == 1:
                # Pass back as dict
                self.memoize[k] = dict(zip(self._data.dtype.names, self._data[mt][0]))
                # return self._data[mt][0]  # Pass back as a Row not a Table
            elif isinstance(k, tuple):
                self.memoize[k] = self._data[mt]
            else:
                raise ValueError(
                    '{:s}: Multiple lines in the list'.format(self.__class__))
            # Finish
            tmp = self.memoize[k].copy()
        return tmp

    # Printing
    def __repr__(self):
        if len(self.sort_by) == 1:
            sort_by = self.sort_by[0]
        else:
            sort_by = self.sort_by
        return '<LineList: {:s}; {} transitions sorted by {:s}>'.format(self.list, len(self._data), sort_by)



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

    sort_by : str or list of str, optional
        Keys to sort the underlying data table by. Default is 'wrest'

    redo_extra_columns : bool, optional
        Whether to recalculate extra columns for log(w*f), abundance, and ion_correction.
        Setting this to True is useful if a different abundance, or ionizatation_correction
        is used. Default is False.

    """

    # Init
    def __init__(self, llst_key, verbose=False, closest=False, set_lines=True,
                 use_ISM_table=True, use_cache=True, sort_by='wrest', redo_extra_columns=False):

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
        else:
            self._data = None

        # Memoize
        self.memoize = {}  # To speed up multiple calls

        # Sort
        self.sort_by = sort_by
        #if (self._data is not None) and (sort_by is not None):
        if self._data is not None:
            # redo extra columns?
            self.set_extra_columns_to_datatable(redo=redo_extra_columns)
            # sort the LineList
            self.sortdata(sort_by)

    @property
    def name(self):
        """ Return the transition names
        """
        return np.array(self._data['name'])

    @property
    def wrest(self):
        """ Return the rest wavelengths
        """
        return self._data['wrest']

    @property
    def Z(self):
        """ Return the Z of the transitions
        """
        return self._data['Z']

    @property
    def ion(self):
        """ Return the ionization state of the transitions
        """
        return self._data['ion']

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
                'LineList: load_data: Not ready for this LineList name: {:s}'.format(self.list))

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
                            mt = np.abs(row['wrest'] - wrest) < tol
                            if mt.sum() == 0:
                                newi.append(jj)
                        # Append new ones as Tables (can't stack QTables yet)
                        full_table = vstack([Table(full_table), Table(table[newi])])
                        full_table = QTable(full_table)
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

    def set_extra_columns_to_datatable(self, abundance_type='solar', ion_correction='none',
                                       redo=False):
        """Sets new convenient columns to the self._data QTable. These will be useful
        for sorting the underlying data table in convenient ways, e.g. by expected
        relative strength, abundance, etc.

        * For atomic transitions these new columns include:
            - `ion_name` : HI, CIII, CIV, etc
            - `log(w*f)` : np.log10(wrest * fosc)  # in np.log10(AA)
            - `abundance` : either [`none`, `solar`]
            - `ion_correction` : [`none`]
            - `rel_strength` : log(w*f) + abundance + ion_correction

        * For molecules a different approach is used:
            - `ion_name` : B0-0P, C6-0, etc.
            - `rel_strength`: We have only three arbitrary levels: [1, 50, 100]
                    100 is for Jk={0,1}
                    50 is for Jk={2,3}
                    1 otherwise

        Parameters
        ----------
        abundance_type : str, optional
            Abundance type. Options are:
                'solar' : Use Solar Abundance (in log10) as given
                          by Asplund 2009. (Default)
                'none' : No abundance given, so this column will
                         be filled with zeros.
        ion_correction: str, optional
            Ionization correction. Options are:
                'none' : No correction applied, so this column will
                         be filled with zeros. (Default)
        redo : bool, optional
            Remake the extra columns

        Note
        ----
        This method is only implemented for the following lists: HI, ISM, EUV, Strong.
        Partially implemented for: H2.

        """
        # Avoid redo (especially for caching)
        if ('ion_name' in self._data.keys()) and (not redo):
            return

        if self.list not in ['HI', 'ISM', 'EUV', 'Strong', 'H2']:
            warnings.warn('Not implemented: will not set relative strength for LineList: {}.'.format(self.list))
            return

        # Set ion_name column
        if self.list in ['HI', 'ISM', 'EUV', 'Strong']:
            ion_name = [name.split(' ')[0] for name in self.name]  # valid for atomic transitions
        elif self.list in ['H2']:
            ion_name = [name.split('(')[0] for name in self.name]  # valid for H2
        self._data['ion_name'] = ion_name

        if self.list in ['H2']:
            # we want Jk to be 1 or 0 first
            cond = (self._data['Jk'] <= 1)
            rel_strength = np.where(cond, 100, 1)
            # second level: Jk = {2,3}
            cond = (self._data['Jk'] > 1) & (self._data['Jk'] <= 3)
            rel_strength = np.where(cond, 50, rel_strength)
            self._data['rel_strength'] = rel_strength
            return

        # Set QM strength as MaskedColumn (self._data['f'] is MaskedColumn)
        qm_strength = self._data['f'] * (self._data['wrest'].to('AA').value)
        qm_strength.name = 'qm_strength'
        self._data['log(w*f)'] = np.log10(qm_strength)
        # mask out potential nans
        cond = np.isnan(self._data['log(w*f)'])
        self._data['log(w*f)'].mask = np.where(cond, True, self._data['log(w*f)'].mask)

        # Set Abundance
        if abundance_type not in ['none', 'solar']:
            raise ValueError('set_extra_columns_to_datatable: `abundance type` '
                             'has to be either: `none` or `solar`')
        if abundance_type == 'none':
            self._data['abundance'] = np.zeros(len(self._data))
        elif abundance_type == 'solar':
            from linetools.abund.solar import SolarAbund
            solar = SolarAbund()
            # abund will be masked array as default (all masked out) in
            # case an element is not in SolarAbund()
            abund = np.ma.masked_array(np.zeros(len(self._data)), mask=True)
            for ii in range(len(abund)):
                ion_name = self._data[ii]['ion_name']
                ion_Z = self._data[ii]['Z']
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
        if ion_correction not in ['none']:
            raise ValueError('set_extra_columns_to_datatable: `ion_correction` '
                             'has to be `none`.')
        if ion_correction == 'none':
            self._data['ion_correction'] = np.zeros(len(self._data))  # is in log10 scale, so 0 means no-correction

        # Set relative strength in log10 scale
        self._data['rel_strength'] = self._data['log(w*f)'] + self._data['abundance'] + self._data['ion_correction']

    def sortdata(self, keys, reverse=False):
        """Sort the LineList according to a given key or keys.

        Parameters
        ----------
        keys : str or list of str
            The main key(s) to sort the LineList by
            (e.g. 'wrest' or ['Z', 'log(w*f)']).
        reverse : bool, optional
            If True, the sorting is reversed
            Default is False.

        Note: this is a wrapper to astropy.table.table.sort()
        """
        # define the sorting key(s) as list
        if isinstance(keys, (str, basestring)):
            keys = [keys]

        # if key is 'as_given', leave it as is
        if keys[0] == 'as_given':
            return

        # sort
        self._data.sort(keys)

        # reverse?
        if reverse:
            self._data.reverse()
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
            initialization (i.e. the default list). This is useful for
            changing subsets of lines without the need to initialize a
            different LineList() object. Default is False.
        sort_by : list of str
            Key(s) to sort the lines by. If sort_by == 'as_given', it will
            preserve the sorting as given by `subset`.

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

        # Sort
        tmp = self._data[indices]
        # Return LineList object
        new = LineList(self.list, closest=self.closest, set_lines=False, verbose=self.verbose)
        new._data = tmp
        if sort_by == ['as_given'] or sort_by == 'as_given':
            pass
        else:
            new.sortdata(sort_by)

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
        """ Get all the transitions corresponding to an ion species
        as given by `line`. In other words, for a given single
        line transition, this function returns all transitions from
        the LineList containing that line.

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

        if self.list not in ['HI', 'ISM', 'EUV', 'Strong', 'H2']:
            warnings.warn('Not implemented for LineList: {}.'.format(self.list))
            return

        if isinstance(line, (str, basestring)):  # Name
            line = line.split(' ')[0]  # keep only the first part of input name
        elif isinstance(line, Quantity):  # Rest wavelength (with units)
            data = self.__getitem__(line)
            return self.all_transitions(data['name'])
        else:
            raise ValueError('Not prepared for this type.')

        if line == 'unknown':
            return self.unknown_line()
        if self.list in ['H2']:

            cond = self._data['ion_name'] == line.split('(')[0]
            tbl = self._data[cond]  # cond is a masked boolean array
            if len(tbl) > 1:
                return tbl
            else:  # this should be always len(tbl)==1 because line was found
                return self.from_qtable_to_dict(tbl)

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
                #cond = tbl['Ej'] == Ej
                cond = np.array([name1.split(' ')[0] == line for name1 in tbl['name']])
                tbl = tbl[cond]
                tbl.sort(['Ej','wrest'])
                # For hydrogen/deuterium this contains deuterium/hydrogen;
                # so let's get rid of them
                #if (line == 'HI') or (line == 'DI'):
                #    names = np.array(tbl['name'])
                #    cond = np.array([l.startswith(line) for l in names])
                #    tbl = tbl[cond]
                if len(tbl) > 1:
                    return tbl
                else:  # this should be always len(tbl)==1 because Z is not None
                    return self.from_qtable_to_dict(tbl)
            else:
                raise ValueError(
                    'Line {} not found in the LineList: {}'.format(line, self.list))

    def strongest_transitions(self, line, wvlims, n_max=3, verbose=False):
        """ Find the strongest transition for an ion

        For a given single line transition, this function returns
        the n_max strongest transitions of the ion species found in
        the LineList, within the wavelength range wlims.

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
        good_linelists = ['HI', 'ISM', 'EUV', 'Strong']
        if self.list not in good_linelists:
            warnings.warn('Not implemented for LineList: {}.'.format(self.list))
            return

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
            return self.from_qtable_to_dict(data[cond])
        else:
            # remove transitions out of range
            data = data[cond]
            # sort by relative strength
            sorted_inds = np.argsort(data['rel_strength'])
            # reverse sorted indices, so strongest get first
            sorted_inds = sorted_inds[::-1]
            # sort using sorted_inds
            data = data[sorted_inds]
            # keep only the first n_max or less
            if n_max is not None:
                data = data[:n_max]
            if len(data) == 1:  # Only 1 case from a QTable format; return a dictionary
                return self.from_qtable_to_dict(data)
            else:
                return data

    def available_transitions(self, wvlims, n_max_tuple=None, min_strength=0.):
        """ Find the strongest transitions in a wavelength interval.

        For a given wavelength range, wvlims = (wv_min, wv_max), this
        function retrieves the n_max_tuple strongest transitions per
        each ion species in the LineList available at such a
        wavelength range and having strength larger than min_strength.
        Strength is defined as log10(wrest * fosc * abundance * ion_correction).
        The output is sorted by strength of the strongest available transition
        of the ion species.

        Parameters
        ----------
        wvlims : tuple of Quantity
            Wavelength range, e.g. wvlims = (1100 * u.AA, 3200 * u.AA)
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

        if self.list not in ['HI', 'ISM', 'EUV', 'Strong']:
            warnings.warn('Not implemented for LineList: {}.'.format(self.list))
            return

        if not all((isinstance(n, int) or (n is None)) for n in [n_max_tuple]):
            raise SyntaxError(
                'n_max_tuple must be integer when given!')
        if isinstance(min_strength, (float,int)):
            pass
        else:
            raise SyntaxError('min_strength must be a float value')

        # Identify unique ion_names (e.g. HI, CIV, CIII)
        unique_ion_names = np.unique(self._data['ion_name'])

        # obtain the strongest available transition of a given unique ion species
        transition_name = []
        strength = []
        for ion in unique_ion_names:  # This loop is necessary to have a non trivial but convenient order in the final output
            aux = self.strongest_transitions(ion, wvlims, n_max=1)  # only the strongest available

            if aux is not None:
                assert isinstance(aux, dict)  # this should always be True given n_max=1
                transition_name += [aux['name']]
                strength += [aux['rel_strength']]
        if len(transition_name) == 0:
            # no matches
            return None

        # create auxiliary Table
        unique = Table()
        unique['name'] = transition_name
        unique['rel_strength'] = strength

        # get rid of those below the min_strength threshold
        cond = unique['rel_strength'] >= min_strength
        unique = unique[cond]
        if len(unique) < 1:
            return None

        # sort by strength
        unique.sort(['rel_strength'])
        unique.reverse()
        # Table unique is now sorted by strength, with only
        # 1 entry per ion species

        # Create output data table adding up to n_max_tuple per ion species
        output = Table()
        for i, row in enumerate(unique):
            name = row['name']
            aux = self.strongest_transitions(name, wvlims, n_max=n_max_tuple)

            # need to deal with dict vs QTable format now
            if isinstance(aux, dict):
                aux = self.from_dict_to_qtable(aux)

            # convert to Table because QTable does not like vstack
            output = vstack([output, Table(aux)])

        # Deal with output formatting now
        # if len==1 return dict
        if len(output) == 1:
            return self.from_qtable_to_dict(output)
        else:
            return QTable(output)

    def from_dict_to_qtable(self, a):
        """Converts dictionary `a` to its QTable version.

        Parameters
        ----------
        a : dict
            The input dictionary to be converted to a
            QTable. The resulting QTable will have the
            same keys as self._data

        Returns
        -------
        A QTable of 1 Row, with filled with the data from
        the input dictionary.

        """

        if isinstance(a, dict):
            pass
        else:
            raise ValueError('Input has to be a dictionary.')

        keys = self._data.keys()

        # Create a QTable with same shape as self._data
        tab = QTable(data=self._data[0])
        # re-write the value elements
        for key in keys:
            tab[0][key] = a[key]
        return tab

    def from_qtable_to_dict(self, tab):
        """Converts QTable `tab` to its dictionary version.
        An error is raised if len(tab) > 1.

        Parameters
        ----------
        tab : QTable or Table
            The table to be converted to a dictionary.
            It has to be of length == 1!

        Returns
        -------
        A dictionary with the same keys and values of the
        input table.

        """

        if not isinstance(tab, (QTable, Table)):
            raise ValueError('Input has to be QTable or Table.')
        elif len(tab) != 1:
            raise ValueError('Input has to be of len(tab) == 1.')

        a_dict = dict()
        for key in tab.keys():
            a_dict[key] = tab[0][key]
        return a_dict

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
        if len(self._data) > 1:
            s = '<LineList: {:s}; {} transitions sorted by {}.>'.format(self.list, len(self._data), self.sort_by)
        else:
            s = '<LineList: {:s}; {} transition.>'.format(self.list, len(self._data))
        return s



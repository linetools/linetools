""" Contains the LineList class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np

from pkg_resources import resource_filename

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.table import Table, hstack
import warnings
import pdb

import imp

lt_path = imp.find_module('linetools')[1]

# from xastropy.xutils import xdebug as xdb

CACHE = {'full_table': {}, 'data': {}}

from linetools.lists import parse as lilp
from linetools.lists import utils as lilu

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
        * 'EUV'     :: EUV lines (for CASBAH project);  Limited X-ray lines too
        * 'Galaxy'  :: Lines typically identified in galaxy spectra
        * 'AGN'     :: Key AGN lines

    verbose : bool, optional
        Give info galore if True

    sort_by : str or list of str, optional
        Keys to sort the underlying data table by. Default is 'wrest'

    redo_extra: bool, optional
        Whether to recalculate extra columns for log(w*f), abundance, and ion_correction.
        Setting this to True is useful if a different abundance, or ionizatation_correction
        is used. Default is False.

    """

    # Init
    def __init__(self, llst_key, verbose=False, closest=False, set_lines=True,
                 use_cache=True, sort_by='wrest', redo_extra=False):

        # Error catching
        if not isinstance(llst_key, basestring):
            raise TypeError('LineList__init__: Wrong type for LineList input')

        # Save
        self.list = llst_key

        # Take closest line?
        self.closest = closest
        self.verbose = verbose

        # Load Data
        self.load_data()
        '''
        if not use_ISM_table or llst_key not in ('ISM', 'HI', 'Strong'):
            self.load_data(use_cache=use_cache)
        '''

        if set_lines:
            # Set lines for use (from defined LineList)
            # This sets self._data
            self.set_lines(verbose=verbose, use_cache=use_cache)
        else:
            self._data = None

        # Memoize
        self.memoize = {}  # To speed up multiple calls

        # Sort
        self.sort_by = sort_by
        # Extras
        self._extra_table = Table(masked=True)
        if self._data is not None:
            # redo extra columns?
            self._extra_table['Id_ex'] = self._data['Id']
            self.make_extra_table(redo=redo_extra)
            # sort the LineList
            self.sortdata(sort_by)

    @property
    def name(self):
        """ Return the transition names
        """
        return np.array(self._data['name'])

    @property
    def wrest(self):
        """ Return the rest wavelengths as a Quantity array
        """
        return Quantity(self._data['wrest'])

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

    def load_data(self):
        data_file = resource_filename('linetools', 'data/lines/linelist.ascii')
        # Read
        self._fulltable = Table.read(data_file, format='ascii.ecsv')

    '''
    def load_data(self, use_ISM_table=True, tol=1e-3, use_cache=True):
        """Grab the data for the lines of interest
        Also load into CACHE

        Parameters
        ----------
        use_ISM_table : bool, optional
        tol : float, optional
          Tolerance for matching wavelength in AA
        use_cache : bool, optional

        Returns
        -------

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
            print('LineList: load_data: Not ready for this LineList name: {:s}'.format(self.list))
            raise ValueError('LineList:  Available options are --  ISM, Strong, EUV, H2, CO, HI, Galaxy')
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
                        for jj, row in enumerate(table):
                            mt = np.abs(row['wrest'] - wrest) < tol
                            if mt.sum() == 0:
                                newi.append(jj)
                        # Append
                        try:
                            full_table = vstack([full_table, table[newi]])
                        except NotImplementedError:
                            pdb.set_trace()
                    # Save to avoid repeating
                    all_func.append(func)

        # Save
        self._fulltable = full_table.copy()

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
    '''

    #####
    def set_lines(self, verbose=True, use_cache=True):
        """ Parse the lines of interest

        Parameters
        ----------
        verbose : bool, optional
        use_cache : bool, optional
          cache the linelist for faster repeat performance
        """
        global CACHE
        key = self.list
        if use_cache and key in CACHE['data']:
            self._data = CACHE['data'][key]
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
        elif self.list == 'AGN':
            set_flags.append('fAGN')
        else:
            raise ValueError(
                'set_lines: Not ready for this: {:s}'.format(self.list))

        # Deal with Defined sets
        if len(set_flags) > 0:
            # Read standard file
            set_data = lilp.read_sets()
            # Speed up
            wrest = self._fulltable['wrest']  # Assuming Angstroms
            for sflag in set_flags:
                gdset = np.where(set_data[sflag] == 1)[0]
                # Match to wavelengths
                for igd in gdset:
                    mt = np.where(np.abs(set_data['wrest'][igd] - wrest) < 9e-5)[0]
                    if len(mt) == 1:
                        self._fulltable['name'][mt[0]] = set_data['name'][igd]
                        indices.append(mt[0])
                    elif len(mt) > 1:
                        #
                        #wmsg = 'WARNING: Multiple lines with wrest={:g}'.format(
                        #    set_data['wrest'][igd])
                        #warnings.warn(wmsg)
                        #warnings.warn(
                        #    'Taking the first entry. Maybe use higher precision.')
                        self._fulltable['name'][mt[0]] = set_data['name'][igd]
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
        tmp_tab['Id'] = np.arange(len(tmp_tab)).astype(int)

        # Finish
        self._data = tmp_tab
        CACHE['data'][key] = self._data


    def make_extra_table(self, abundance_type='solar', ion_correction='none',
                                       redo=False):
        """Build an additional table that is parallel to self._data that
        includes convenient columns. These will be useful
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
        if ('ion_name' in self._extra_table.keys()) and (redo is False):
            return

        if self.list not in ['HI', 'ISM', 'EUV', 'Strong', 'H2']:
            warnings.warn('Not implemented: will not set relative strength for LineList: {}.'.format(self.list))
            return

        # Set ion_name column
        ion_name = np.array([' '*20]*len(self.name)).astype(str)
        if self.list in ['HI', 'ISM', 'EUV', 'Strong']:
            for kk,name in enumerate(self.name):  # valid for atomic transitions
                ion_name[kk] = name.split(' ')[0]
        elif self.list in ['H2']:
            for kk,name in enumerate(self.name):  # H2
                ion_name[kk] = name.split('(')[0]
        self._extra_table['ion_name'] = ion_name

        if self.list in ['H2']:
            # we want Jk to be 1 or 0 first
            cond = (self._data['Jk'] <= 1)
            rel_strength = np.where(cond, 100, 1)
            # second level: Jk = {2,3}
            cond = (self._data['Jk'] > 1) & (self._data['Jk'] <= 3)
            rel_strength = np.where(cond, 50, rel_strength)
            self._extra_table['rel_strength'] = rel_strength
            return

        # Set QM strength as MaskedColumn (self._data['f'] is MaskedColumn)
        qm_strength = self._data['f'] * self._data['wrest']
        qm_strength.name = 'qm_strength'
        self._extra_table['log(w*f)'] = np.zeros_like(qm_strength)
        gdqm = qm_strength > 0.
        self._extra_table['log(w*f)'][gdqm] = np.log10(qm_strength[gdqm])
        # mask out potential nans
        self._extra_table['log(w*f)'].mask = ~gdqm

        # Set Abundance
        if abundance_type not in ['none', 'solar']:
            raise ValueError('set_extra_columns_to_datatable: `abundance type` has to be either: `none` or `solar`')
        if abundance_type == 'none':
            self._extra_table['abundance'] = np.zeros(len(self._data))
        elif abundance_type == 'solar':
            from linetools.abund.solar import SolarAbund
            solar = SolarAbund()
            abund = np.ma.masked_array(np.zeros(len(self._data)), mask=True)
            for row in solar._data:
                Zmt = self._data['Z'] == row['Z']
                abund[Zmt] = row['Abund']
                abund.mask[Zmt] = False
            # Deuterium is special
            fchar = np.array([ion_name[0] for ion_name in self._extra_table['ion_name'].data])
            DI = fchar == 'D'
            abund[DI] = solar['D']
            # Finish
            self._extra_table['abundance'] = abund

        # Set ionization correction
        if ion_correction not in ['none']:
            raise ValueError('set_extra_columns_to_datatable: `ion_correction` '
                             'has to be `none`.')
        if ion_correction == 'none':
            self._extra_table['ion_correction'] = np.zeros(len(self._data))  # is in log10 scale, so 0 means no-correction

        # Set relative strength in log10 scale
        self._extra_table['rel_strength'] = self._extra_table['log(w*f)'] + self._extra_table['abundance'] + self._extra_table['ion_correction']

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

        # hstack for convenience
        if len(self._extra_table) == len(self._data):
            flg_extra = True
            dkeys = self._data.keys()
            ekeys = self._extra_table.keys()
            dtbl = hstack([self._data, self._extra_table], join_type='exact')
        else:
            flg_extra = False
            dtbl = self._data

        # sort
        dtbl.sort(keys)

        # reverse?
        if reverse:
            dtbl.reverse()
        self.sort_by = keys

        # Finish
        if flg_extra:
            self._data = dtbl[dkeys]
            self._extra_table = dtbl[ekeys]
        else:
            self._data = dtbl

    def subset_lines(self, subset, reset_data=False, verbose=False, sort_by=['wrest']):
        """ Select a user-specific subset of the lines from the LineList
        Code does *not* raise an error if a requested line is not present

        Parameters
        ----------
        subset : list or Quantity array (of Quantity or str)
            List of wrest or names for lines to use (drawn from current
            LineList object) e.g. (['HI 1215', 'CIV 1548'] or
            [1215.67, 1548.195] * u.AA
        reset_data : bool, optional
            Reset self._data Table based on the original list at the
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
        if not isinstance(subset, (list, Quantity)):
            raise ValueError('subset_lines: the input subset must be a list or Quantity array!')

        # Reset _data (useful for changing subsets)
        if reset_data:
            self.set_lines(verbose=False)

        indices = []
        if isinstance(subset, Quantity):  # wrest
            wrest = self._data['wrest'].to('AA').value
            # Generate dummy 2D arrays
            subset2d = np.outer(subset.to('AA').value, np.ones_like(wrest))
            wrest2d = np.outer(np.ones(len(subset)), wrest)
            diff = np.abs(subset2d-wrest2d)
            indices = np.where(diff < 1e-4)[1].tolist()
        elif isinstance(subset[0], (basestring)):  # Names
            # Using sets
            names = set(self._data['name'])
            inter = set(subset).intersection(names)  # But these aren't ordered the same
            # For the indices
            ind_dict = dict((k,i) for i,k in enumerate(self._data['name']))
            indices = [ind_dict[x] for x in inter]
        else:
            raise ValueError('Not ready for this `subset` type yet.')

        # Sort
        tmp = self._data[indices]
        # Return LineList object
        new = LineList(self.list, closest=self.closest,
                       set_lines=False, verbose=self.verbose)
        new._data = tmp
        if sort_by == ['as_given'] or sort_by == 'as_given':
            if isinstance(subset, Quantity):  # wrest
                pass
            elif isinstance(subset[0], (basestring)):  # Names
                isort = []
                names = list(new._data['name'])
                for name in subset:
                    isort.append(names.index(name))
                new._data = new._data[isort]
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

    def all_transitions(self, line, debug=False):
        """ Get all the transitions corresponding to an ion species
        as given by `line`. In other words, for a given single
        line transition, this function returns all transitions from
        the LineList containing that line.

        Parameters
        ----------
        line : str or Quantity or tuple
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII',
            1215.6700*u.AA). When string contains spaces it only
            considers the first part of it, so 'HI' and 'HI 1215' and
            'HI 1025' are all equivalent. To retrieve an unknown line
            use string 'unknown'.

            If tuple -- Zion, e.g. (6,2)

        Returns
        -------
        dict or Table
            dict if only 1 transition found, otherwise Table.
            If dict, it will have units!

        """

        if self.list not in ['HI', 'ISM', 'EUV', 'Strong', 'H2']:
            warnings.warn('Not implemented for LineList: {}.'.format(self.list))
            return

        if isinstance(line, (str, basestring)):  # Name
            line = line.split(' ')[0]  # keep only the first part of input name
        elif isinstance(line, (Quantity, tuple)):  # Rest wavelength (with units)
            data = self.__getitem__(line)
            if isinstance(line, tuple):
                data=data[0]
            return self.all_transitions(data['name'])
        else:
            raise ValueError('Not prepared for this type.')

        if line == 'unknown':
            return self.unknown_line()
        if self.list in ['H2']:
            cond = self._extra_table['ion_name'] == line.split('(')[0]
            tbl = self._data[cond]  # cond is a masked boolean array
            if len(tbl) > 1:
                return tbl
            else:  # this should be always len(tbl)==1 because line was found
                return lilu.from_table_to_dict(tbl)

        else:
            Z = None
            indices = np.where(self._extra_table['ion_name'] == line)[0]
            if len(indices) == 0:
                raise ValueError(
                    'Line {} not found in the LineList: {}'.format(line, self.list))
            else:
                idx = indices[0]
            Z = self._data['Z'][idx]  # atomic number
            ie = self._data['ion'][idx]  # ionization state
            Ej = self._data['Ej'][idx]  # Energy of lower level
            tbl = self.__getitem__((Z, ie))
            # Make sure the lower energy level is the same too
            if debug:
                pdb.set_trace()
            cond = tbl['Ej'] == Ej
            tbl = tbl[cond]
            tbl.sort(['Ej','wrest'])
            # For hydrogen/deuterium this contains deuterium/hydrogen;
            # so let's get rid of them
            if (line == 'HI') or (line == 'DI'):
                names = np.array(tbl['name'])
                cond = np.array([l.startswith(line) for l in names])
                tbl = tbl[cond]
            if len(tbl) > 1:
                return tbl
            else:  # this should be always len(tbl)==1 because Z is not None
                return lilu.from_table_to_dict(tbl)

    def strongest_transitions(self, line, wvlims, n_max=3, verbose=False,
                              debug=False):
        """ Find the strongest transition for a single ion
        available_transitions() finds those for all ions in a wavelength interval

        For a given single line transition, this function returns
        the n_max strongest transitions of the ion species found in
        the LineList, within the wavelength range wlims.

        Parameters
        ----------
        line : str or Quantity or tuple
            Name of line. (e.g. 'HI 1215', 'HI', 'CIII', 'SiII',
            1215.6700*u.AA). When string contains spaces it only
            considers the first part of it, so 'HI' and 'HI 1215' and
            'HI 1025' are all equivalent. To retrieve an unknown line
            use string 'unknown'.

            If tuple -- Zion, e.g. (6,2)

        wvlims : tuple of Quantity, or Quantity tuple
            Wavelength range, e.g. wvlims = (1100 * u.AA, 3200 * u.AA) or
            wvlims = (1100, 3200) * u.AA
        n_max : int or None, optional
            Maximum number of transitions to retrieve; if n_max = None
            it retrieves all of them

        Returns
        -------
        None (if no transitions are found), dict (if only 1 transition
        found), or Table (if > 1 transitions are found)

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

        # Grab all_transitions related to the line
        data = self.all_transitions(line)
        # condition to be within wvrange
        cond = (Quantity(data['wrest']) >= wvlims[0]) & (Quantity(data['wrest']) <= wvlims[1])
        if np.sum(cond) == 0:
            if verbose:
                print(
                    '[strongest_transitions] Warning: no transitions found within wvlims; returning None')
            return None
        elif isinstance(data, dict):  # Only 1 case from a dict format
            return data
        elif np.sum(cond) == 1:  # only 1 case from a Table format
            return lilu.from_table_to_dict(data[cond])
        else:
            # remove transitions out of range
            data = data[cond]
            # sort by relative strength
            idx = data['Id'].data
            sorted_inds = np.argsort(self._extra_table['rel_strength'][idx])
            # reverse sorted indices, so strongest get first
            sorted_inds = sorted_inds[::-1]
            # sort using sorted_inds
            data = data[sorted_inds]
            # keep only the first n_max or less
            if n_max is not None:
                data = data[:n_max]
            if len(data) == 1:  # Only 1 case from a Table format; return a dictionary
                return lilu.from_table_to_dict(data)
            else:
                return data

    def available_transitions(self, wvlims, n_max_tuple=None, min_strength=0.):
        """ Find the strongest available transitions in a wavelength interval.

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
        dict (if only 1 transition found) or Table (if > 1
        transitions are found) or None (if no transition is found)
        This is an hstack of self._data and self._extra_table
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

        # Cut on wavelength first!
        gdwv = (self.wrest>wvlims[0]) & (self.wrest < wvlims[1])
        if np.any(gdwv):
            tmp_data = self._data[gdwv]
            extras = self._extra_table[gdwv]
            data = hstack([tmp_data,extras])
        else:
            return None

        # Sort
        data.sort(['ion_name','rel_strength'])
        data.reverse()


        # Identify unique ion_names (e.g. HI, CIV, CIII)
        # This is WRONG for FeII*
        unique_ion_names, idx = np.unique(data['ion_name'], return_index=True)
        unique= data[idx]

        '''
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
        pdb.set_trace()
        unique = Table()
        unique['name'] = transition_name
        unique['rel_strength'] = strength
        '''

        # get rid of those below the min_strength threshold
        cond = unique['rel_strength'] >= min_strength
        unique = unique[np.where(cond == True)[0]] # Deals with masking
        if len(unique) < 1:
            return None

        # sort by strength
        unique.sort(['rel_strength'])
        unique.reverse()
        # Table unique is now sorted by strength, with only
        # 1 entry per ion species

        # Create output data table adding up to n_max_tuple per ion species
        indices = []
        for ion_name in unique['ion_name']:
            mt = np.where(data['ion_name'] == ion_name)[0]
            indices += mt[0:n_max_tuple].tolist()
        output = data[np.array(indices)]

        # Deal with output formatting now
        # if len==1 return dict
        if len(output) == 1:
            return lilu.from_table_to_dict(output)
        else:
            return output


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
          Table (tuple when more than 1 lines are found)
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
                tmp2 = {}
                for name in self._data.dtype.names:  # Captures units
                    if self._data[name].unit is None:
                        tmp2[name] = self._data[name][mt][0]
                    else:
                        tmp2[name] = self._data[name][mt][0] * self._data[name].unit
                self.memoize[k] = tmp2.copy()
                # return self._data[mt][0]  # Pass back as a Row not a Table
            elif isinstance(k, tuple):
                self.memoize[k] = self._data[mt]
            else:
                raise ValueError(
                    '{:s}: Multiple lines in the list with your input.  Give a more unique input or change the tol.'.format(self.__class__.__name__))
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



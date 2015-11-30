""" Class for and AbslineSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import json
import copy
from abc import ABCMeta
import pdb

from astropy.io import ascii 
from astropy import units as u
from astropy.table import QTable, Column
from astropy.units.quantity import Quantity

from linetools.spectra import io as lsio
from linetools.isgm.abssystem import AbsSystem


class AbslineSurvey(object):
    """
    Class for a survey of absorption line systems.

    Attributes
    ----------
        abs_type : str, unicode
          Type of Absorption system (DLA, LLS)
        ref : str
          Reference(s) to the Survey
    """

    __metaclass__ = ABCMeta

    @classmethod
    def from_flist(cls, flist, tree=None, **kwargs):
        """ Read from list of .dat files (historical JXP format)

        Parameters
        ----------
        flist : str
          ASCII file including list of .dat files
        tree : str, optional
          Path to .dat files
        kwargs :
          Passed to __init__
        """
        if tree is None:
            tree = ''
        # Load up (if possible)
        data = ascii.read(tree+flist, data_start=0,
                          guess=False, format='no_header')
        slf = cls(**kwargs)
        slf.tree = tree
        slf.flist = flist

        # Load up
        slf.dat_files = list(data['col1'])
        # Generate AbsSys list
        for dat_file in slf.dat_files:
            slf._abs_sys.append(set_absclass(slf.abs_type).from_datfile(dat_file, tree=slf.tree))
        print('Read {:d} files from {:s} in the tree {:s}'.format(
            slf.nsys, slf.flist, slf.tree))

        return slf

    @classmethod
    def from_sfits(cls, summ_fits, **kwargs):
        """Generate the Survey from a summary FITS file

        Handles SPEC_FILES too.

        Parameters
        ----------
        summ_fits : str
          Summary FITS file
        **kwargs : dict
          passed to __init__
        """
        # Init
        slf = cls(**kwargs)
        # Read
        systems = QTable.read(summ_fits)
        nsys = len(systems)
        # Dict
        kdict = dict(NHI=['NHI', 'logNHI'],
                     sig_NHI=['sig(logNHI)'],
                     name=['Name'], vlim=['vlim'],
                     zabs=['Z_LLS'], zem=['Z_QSO'],
                     RA=['RA'], Dec=['DEC', 'Dec'])
        # Parse the Table
        inputs = {}
        for key in kdict.keys():
            vals, tag = lsio.get_table_column(kdict[key], [systems],idx=0)
            if vals is not None:
                inputs[key] = vals
        # vlim
        if 'vlim' not in inputs.keys():
            default_vlim = [-500, 500.]* u.km / u.s
            inputs['vlim'] = [default_vlim]*nsys
        # Generate
        for kk in range(nsys):
            # Generate keywords
            kwargs = {}
            args = {}
            for key in inputs.keys():
                if key in ['vlim', 'zabs', 'RA', 'Dec']:
                    args[key] = inputs[key][kk]
                else:
                    kwargs[key] = inputs[key][kk]
            # Instantiate
            abssys = set_absclass(slf.abs_type)((args['RA'], args['Dec']), args['zabs'], args['vlim'], **kwargs)
            # spec_files
            try:
                abssys.spec_files += systems[kk]['SPEC_FILES'].tolist()
            except (KeyError, AttributeError):
                pass
            slf._abs_sys.append(abssys)
        # Return
        return slf

    def __init__(self, abs_type, ref=''):
        # Expecting a list of files describing the absorption systems
        """  Initiator

        Parameters
        ----------
        abs_type : str, unicode
          Type of AbsSystem in the Survey, e.g.  MgII, DLA, LLS
        ref : string, optional
          Reference(s) for the survey
        """
        self.abs_type = abs_type
        self.ref = ref
        self._abs_sys = []
        self.mask = None

        # Mask
        self.init_mask()

        # Init
        self.flist = None

    @property
    def nsys(self):
        """ Number of systems
        Returns
        -------
        nsys : int
        """
        return len(self._abs_sys)

    def init_mask(self):
        """ Initialize the mask for abs_sys
        """
        if self.nsys > 0:
            self.mask = np.array([True]*self.nsys)

    def abs_sys(self):
        # Recast as an array
        return lst_to_array(self._abs_sys, mask=self.mask)

    def add_abs_sys(self, abs_sys):
        """ Add an AbsSys to the Survey

        Enables one to add checks

        Parameters
        ----------
        abs_sys : AbsSystem
        """
        assert self.chk_abs_sys(abs_sys)
        # Might check to see if a duplicate exists..

        # Append
        self._abs_sys.append(abs_sys)

    def chk_abs_sys(self, abs_sys):
        """ Preform checks on input abs_sys

        Parameters
        ----------
        abs_sys : AbsSystem

        Returns
        -------
        bool

        """
        if not isinstance(abs_sys, AbsSystem):
            raise IOError("Must be an AbsSystem object")
        return True

    def __getattr__(self, k):
        """ Generate an array of attribute 'k' from the AbsSystems

        Mask is applied

        Parameters
        ----------
        k : str
          Attribute

        Returns
        -------
        numpy array
        """
        try:
            lst = [getattr(abs_sys, k) for abs_sys in self._abs_sys]
        except ValueError:
            raise ValueError
        # Recast as an array
        return lst_to_array(lst, mask=self.mask)

    def fill_ions(self, jfile=None):  # This may be overloaded!
        """ Loop on systems to fill in ions

        Parameters
        ----------
        jfile : str, optional
          JSON file containing the information
        """
        if jfile is not None:
            # Load
            with open(jfile) as data_file:    
                ions_dict = json.load(data_file)
            # Loop on systems
            for abs_sys in self._abs_sys:
                abs_sys.get_ions(idict=ions_dict[abs_sys.name])
        else:
            for abs_sys in self._abs_sys:
                # Line list
                if (abs_sys.linelist is None) & (self.linelist is not None):
                    abs_sys.linelist = self.linelist
                #
                abs_sys.get_ions()

    # Get ions
    def ions(self, iZion, skip_null=False):
        """
        Generate a Table of columns and so on
        Restrict to those systems where flg_clm > 0

        Parameters
        ----------
        iZion : tuple
           Z, ion   e.g. (6,4) for CIV
        skip_null : boolean (False)
           Skip systems without an entry, else pad with zeros 

        Returns
        -------
        Table of values for the Survey
        """
        keys = [u'name', ] + self.abs_sys()[0]._ionclms.keys()
        t = copy.deepcopy(self.abs_sys()[0]._ionclms[0:1])
        t.add_column(Column(['dum'], name='name', dtype='<U32'))
        t = t[keys]

        # Loop on systems (Masked)
        for abs_sys in self.abs_sys():
            # Grab
            mt = (abs_sys._ionclms['Z'] == iZion[0]) & (abs_sys._ionclms['ion'] == iZion[1])
            if np.sum(mt) == 1:
                irow = abs_sys._ionclms[mt]
                # Cut on flg_clm
                if irow['flag_N'] > 0:
                    row = [abs_sys.name] + [irow[key] for key in keys[1:]]
                    t.add_row(row)   # This could be slow
                else:
                    if skip_null is False:
                        row = [abs_sys.name] + [0 for key in keys[1:]]
                        t.add_row(row)
            elif np.sum(mt) == 0:
                if skip_null is False:
                    row = [abs_sys.name] + [0 for key in keys[1:]]
                    t.add_row( row )
                continue
            else:
                raise ValueError('Multiple entries...')

        # Return
        return t[1:]

    # Mask
    def update_mask(self, mask, increment=False):
        """ Update the Mask for the abs_sys

        Parameters
        ----------
        mask : array (usually Boolean)
           Mask of systems 
        increment : bool, optional
           Increment the mask (i.e. keep False as False)
        """
        if len(mask) == len(self._abs_sys):  # Boolean mask
            if increment is False:
                self.mask = mask
            else:
                self.mask = self.mask & mask
        else:
            raise ValueError('abs_survey: Needs developing!')

    def __repr__(self):
        if self.flist is not None:
            return '[AbslineSurvey: {:s} {:s}, nsys={:d}, type={:s}, ref={:s}]'.format(
                self.tree, self.flist, self.nsys, self.abs_type, self.ref)
        else:
            return '[AbslineSurvey: nsys={:d}, type={:s}, ref={:s}]'.format(
                self.nsys, self.abs_type, self.ref)


class GenericAbsSurvey(AbslineSurvey):
    """A simple absorption line survey
    """
    def __init__(self, **kwargs):
        AbslineSurvey.__init__(self, 'Generic', **kwargs)


def set_absclass(abstype):
    """Translate abstype into Class

    Parameters
    ----------
    abstype : str
      AbsSystem type, e.g. 'LLS', 'DLA'

    Returns
    -------
    Class name
    """
    from xastropy.igm.abs_sys.dla_utils import DLASystem
    from linetools.isgm.lls import LLSSystem
    from linetools.isgm.abssystem import GenericAbsSystem

    cdict = dict(LLS=LLSSystem, DLA=DLASystem)
    try:
        return cdict[abstype]
    except KeyError:
        return GenericAbsSystem


def lst_to_array(lst, mask=None):
    """ Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Parameters
    ----------
    lst : list
      Should be number or Quantities
    mask : boolean array, optional

    Returns
    -------
    array or Quantity array

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]


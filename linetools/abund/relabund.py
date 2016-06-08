""" Class for handling relative abundances
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str
import numbers

import numpy as np
import pdb

#from astropy.utils.misc import isiterable
from astropy.table import Table

from .solar import SolarAbund
import linetools.abund.elements as ltae

#
class RelAbund(object):
    """Class to handle relative abundances, usually of an AbsSystem

    Parameters
    ----------
    """
    @classmethod
    def from_clm_pair(cls, X, NX, Y, NY, sigNX=0.1, sigNY=0.1,
                      NH=None, **kwargs):
        """ Instantiate from a pair of column density measurements

        Parameters
        ----------
        X : str or int
          Element 1 (e.g. 'C' or 6)
        NX : float
          log10 column of X
        Y : str or int
          Element 2 (e.g. 'C' or 6)
        NY : float
          log10 column of Y
        sigNX : float, optional
          Error in NX
        sigNY : float, optional
          Error in NY
        NH : list, optional
          flag_NH, logNH, sig_logNH

        Returns
        -------
        RelAbund
          Also prints value of [X/Y]

        """
        # Instantiate
        slf = cls(**kwargs)
        # NH
        if NH is None:
            slf.NH = [1,21.0,0.1]
            print('Adopting arbitrary log NH={:f} for [X/H] values'.format(slf.NH[1]))
        else:
            slf.NH = NH
        # Add em in
        slf.add_elmbyclm(1, X, NX, sigNX)
        slf.add_elmbyclm(1, Y, NY, sigNY)
        # Print as we go
        print("Input ratio [{}/{}]={}".format(X,Y,slf[X,Y]))
        return slf

    @classmethod
    def from_ionclm_table(cls, NHI, tbl, low_ions=True, **kwargs):
        """ Generate class from an input table of ionic column densities

        Parameters
        ----------
        NHI : list  [int, float, float]
          flag_NHI, logNHI, sig_logNHI
        tbl : Table
        low_ions : bool, optional
          Generate abundances from low-ions only.  No ionization corrections
          Best for ISM, DLAs only
        kwargs

        Returns
        -------
        slf

        """
        # Checks
        if not low_ions:
            raise IOError("Only coded for low-ions so far")
        tblkeys = tbl.keys()
        for key in ['Z', 'ion', 'flag_N', 'logN', 'sig_logN']:
            if key not in tblkeys:
                raise IOError("Input table must include {:s}".format(key))
        if NHI[0] != 1:
            raise IOError("Not ready for this NHI flag")
        # Start the class
        slf = cls(**kwargs)
        # Store NH
        slf.NHI = NHI
        if low_ions:
            slf.NH = slf.NHI
        # Loop through the input Table
        for row in tbl:
            # Check Ej -- ground-state only
            if 'Ej' in tblkeys:
                if row['Ej'] > 0.:  # Not expecting units
                    continue
            # Skip Hydrogen
            if row['Z'] == 1:
                continue
            # Low-ion?
            elm = slf.elements[row['Z']]
            if row['ion'] == 1:
                if elm.ionenergy[row['ion']-1] < 13.6:
                    continue
            else:
                if (elm.ionenergy[row['ion']-2] > 13.6) or (elm.ionenergy[row['ion']-1] < 13.6):
                    continue
            # Calculate
            XH = row['logN'] - NHI[1] + 12 - slf.solar[row['Z']]
            sigXH = np.sqrt(row['sig_logN']**2 + NHI[2]**2)  # Crude but ok
            # Fill it up
            slf._data[row['Z']] = dict(flag=row['flag_N'], XH=XH, sigXH=sigXH,
                                       sig=row['sig_logN'],  # For relative abundances
                                       )
        # Return
        return slf

    def __init__(self, solar_ref='Asplund2009', verbose=False):
        """
        solar_ref : str, optional
            Reference for the underlying Solar abundances
            'Asplund2009' :: Asplund et al. 2009, ARA&A, 47, 481 meteoritic
            table (several photometric)
        """
        # Init
        self.solar_ref = solar_ref
        self._data = {}  # Nested dict.  Top keys are atomic number (6, 14, 26)
            #  Next dict keys are flag, XH, sigXH, sig

        # Load Solar abundances
        self.solar = SolarAbund(ref=self.solar_ref, verbose=verbose)
        # Load ELEMENTS too
        self.elements = ltae.ELEMENTS

    def add_elmbyclm(self, flag, X, NX, sigNX, clobber=False):
        """ Add an entry with a column density
        Requires NH to be set previously

        Parameters
        ----------
        flag
        X
        NX
        sigNX

        Returns
        -------

        """
        if hasattr(self,'NH') is False:
            raise IOError("This method requires self.NH to be set previously")
        # Setup
        Xint = self.elements[X].number
        if (Xint in self._data.keys()) and (clobber is False):
            print("Not overwriting elm={:d}.  Use clobber to do so".format(Xint))
        XH = NX - self.NH[1] + 12 - self.solar[Xint]
        # Write
        self._data[Xint] = dict(flag=flag, XH=XH,
                                sigXH=np.sqrt(self.NH[2]**2 + sigNX**2),
                                sig=sigNX)

    def table(self, Y=1):
        """ Generate an [X/Y] table from the dict.  Default is [X/H]

        Parameters
        ----------
        Y : int or str, optional
          Relative abundance

        Returns
        -------

        """
        # Init
        if isinstance(Y,basestring):
            Yint = self.elements[Y].number
        elif isinstance(Y,numbers.Integral):
            Yint = Y
        else:
            raise IOError("Bad Y input {}".format(Y))
        Yc = self.elements[Yint].symbol
        #
        clms = ['flag', '[X/{:s}]'.format(Yc), 'sig([X/{:s}])'.format(Yc)]
        dkeys = ['flag', 'val', 'sig']
        lists = [[] for x in xrange(len(clms))]
        # List it
        Zlist, nlist = [], []
        for key in self._data.keys():
            # Skip itself
            if key == Yint:
                continue
            # Z, name
            Zlist.append(key)
            nlist.append(self.elements[key].symbol)
            #
            XYdict = self[key,Yint]
            for jj,dkey in enumerate(dkeys):
                lists[jj].append(XYdict[dkey])
        # Generate the Table
        tbl = Table()
        tbl['Z'] = Zlist
        tbl['Name'] = nlist
        for jj,clm in enumerate(clms):
            tbl[clm] = lists[jj]
        # Sort
        tbl.sort('Z')
        # Meta
        tbl['flag'].meta = {1:'Value', 2:'Lower limit', 3:'Upper Limit'}
        # Return
        return tbl

    def __getitem__(self, k):
        """ Return XY abundance relative to solar as a dict
        given an element or pair of elements.
        If k is a scalar, [X/H] is returned.  Otherwise [X/Y]
 
        Parameters
        ----------
        k : int or str or list/tuple
          * int -- Atomic number (6)
          * str -- Element name (e.g. 'C')
          * tuple/list -- (X,Y) , e.g. (14,26) for [Si/Fe]

        Returns
        -------
        odict_XY : dict
           * 'flag'-- -1=NG, 1=Good value, 2=Lower limit, 3=Upper limit
           * 'val' -- [X/Y]
           * 'sig' -- sigma([X/Y])  rough estimate
        """
        flag_XH = True
        if isinstance(k, (numbers.Integral, basestring)): # XH
            Xint = self.elements[k].number
            XHdict = self._data[Xint]
        elif isinstance(k, (tuple,list)):  # XY
            Xint = self.elements[k[0]].number
            XHdict = self._data[Xint]
            if k[1] not in [1,'H']:
                flag_XH = False
                Yint = self.elements[k[1]].number
                YHdict = self._data[Yint]
                XYdict = dict(flag=0, val=XHdict['XH']-YHdict['XH'],
                              sig=np.sqrt(XHdict['sig']**2 + YHdict['sig']**2))
                # Parse flags
                if YHdict['flag'] == 1:
                    XYdict['flag'] = XHdict['flag']
                elif YHdict['flag'] == 2:  # Lower limit on Y
                    if XHdict['flag'] in [1,3]:
                        XYdict['flag'] = 3
                    else:
                        XYdict['flag'] = -1
                elif YHdict['flag'] == 3:  # Upper limit on Y
                    if XHdict['flag'] in [1,2]:
                        XYdict['flag'] = 2
                    else:
                        XYdict['flag'] = -1
        else:
            raise IndexError('Not prepared for this type of input', k)

        # Return
        if flag_XH:
            return dict(flag=XHdict['flag'], val=XHdict['XH'], sig=XHdict['sigXH'])
        else:
            return XYdict


    def __repr__(self):
        return ('<{:s}:>'.format(self.__class__.__name__))

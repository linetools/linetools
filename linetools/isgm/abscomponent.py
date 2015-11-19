"""
#;+ 
#; NAME:
#; spectralline
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for SpectralLine class
#;   23-Jun-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import pdb
import copy, imp

from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.analysis import utils as lau
from linetools.analysis import absline as laa
from linetools.lists.linelist import LineList

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# Class for Components
class AbsComponent(object):
    """Class for a spectral line.  Emission or absorption 
    Attributes
    ----------
    coord : SkyCoord
    Zion : tuple 
       Atomic number, ion -- (int,int) 
       e.g. (8,1) for OI
    zcomp : float
       Component redshift
    vlim : Quantity array
       Velocity limits of the component
       e.g.  [-300,300]*u.km/u.s
    A : int
        Atomic mass -- used to distinguish isotopes
    Ej : Quantity
       Energy of lower level (1/cm)
    """
    @classmethod
    def from_abslines(cls, abslines):
        """Instantiate from a list of AbsLine objects
        Parameters
        ----------
        abslines : list 
          List of AbsLine objects
        """
        # Check
        if not isinstance(abslines,list):
            raise IOError('Need a list of AbsLine objects')
        # Instantiate with the first line
        init_line = abslines[0]
        slf = cls( init_line.attrib['coord'],
           (init_line.data['Z'],init_line.data['ion']),
           init_line.attrib['z'], init_line.analy['vlim']) 
        slf._abslines.append(init_line)
        # Append with component checking
        if len(abslines) > 1:
            for absline in abslines[1:]:
                slf.add_absline(absline)
        # Return
        return slf

    # Initialize with wavelength
    def __init__(self, radec, Zion, z, vlim, Ej=0./u.cm, A=None, Ntup=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
          (RA,DEC) in deg
        Zion : tuple 
          Atomic number, ion -- (int,int) 
             e.g. (8,1) for OI
        z : float
          Absorption redshift
        vlim : Quantity array
          Velocity limits of the component
          e.g.  [-300,300]*u.km/u.s
        A : int, optional
          Atomic mass -- used to distinguish isotopes
        Ntup : tuple
          (int,float,float)
          (flgN,logN,sigN)
          flgN : Flag describing N measurement
          logN : log10 N column density
          sigN : Error in log10 N
        Ej : Quantity, optional
           Energy of lower level (1/cm)
        """

        # Required
        if isinstance(radec,(tuple)):
            self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        elif isinstance(radec,SkyCoord):
            self.coord = radec
        self.Zion = Zion
        self.zcomp = z
        self.vlim = vlim

        # Optional
        self.A = A
        self.Ej = Ej
        if Ntup is not None:
            self.flgN = Ntup[0]
            self.logN = Ntup[1]
            self.sigN = Ntup[2]

        # Other
        self._abslines = []

    def add_absline(self,absline):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.

        Parameters
        ----------
        absline : AbsLine
        """
        # Perform easy checks
        test = self.coord.separation(absline.attrib['coord']) < 0.1*u.arcsec
        test = test & self.Zion[0] == absline.data['Z']
        test = test & self.Zion[1] == absline.data['ion']
        test = test & bool(self.Ej == absline.data['Ej'])
        # Now redshift/velocity
        # Isotope
        if self.A is not None:
            raise ValueError('Not ready for this yet')
        # Append?
        if test:
            self._abslines.append(absline)
        else:
            print('Input absline with wrest={:g} does not match component rules'.format(absline.wrest))
            print('Not appending')

    # Output
    def __repr__(self):
        return ('[AbsComponent: {:s} {:s}, Zion=({:d},{:d}), z={:g}]'.format(
                self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                self.Zion[0], self.Zion[1], self.zcomp))


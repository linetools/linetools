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

# class SpectralLine(object):
# class AbsLine(SpectralLine):
# class AbsComponens(AbsLine):

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
    A : int, optional
      Atomic mass -- used to distinguish isotopes
    """
    #@classmethod
    #def from_table(cls, table, dispersion_column='dispersion'):

    # Initialize with wavelength
    def __init__(self, radec, Zion, z, vlim, A=None, Ntup=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple
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
        """

        # Required
        self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        self.Zion = Zion
        self.zcomp = z
        self.vlim = vlim

        # Optional
        self.A = A
        if Ntup is not None:
            self.flgN = Ntup[0]
            self.logN = Ntup[1]
            self.sigN = Ntup[2]

    # Output
    def __repr__(self):
        return ('[AbsComponent: {:s} {:s}, Zion=({:d},{:d}), z={:g}]'.format(
                self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                self.Zion[0], self.Zion[1], self.zcomp))


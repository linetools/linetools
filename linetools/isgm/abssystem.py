""" Class for  absorption systems
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import warnings
from abc import ABCMeta

from astropy import units as u
from astropy.table import QTable
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu


class AbsSystem(object):
    """
    Class for an absorption line system

    Parameters
    ----------
    abs_type : str or unicode
      Type of absorption system, e.g. MgII, DLA, LLS
    radec : tuple or coordinate
        RA/Dec of the sightline or astropy.coordinate
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
    zem : float, optional
      Emission redshift of background source
    NHI :  float, optional
      Log10 of the HI column density
    sig_NHI :  np.array(2), optional
      Log10 error of the HI column density (-/+)
    name : str, optional
      Name for the system

    Attributes
    ----------
        abs_type : str
        coord : SkyCoord
            RA/Dec of the sightline
        zabs : float
          Absorption redshift
        zem : float
          Emission redshift of background source
        vlim : Quantity array (2)
          Velocity limits of the system
        NHI :  float
          Log10 of the HI column density
        sig_NHI :  np.array(2)
          Log10 error of the HI column density (-/+)
        MH :  float
          Metallicity (log10)
        name : str
            Name of the system
    """

    __metaclass__ = ABCMeta

    @classmethod
    def from_components(cls, components):
        """Instantiate from a list of AbsComponent objects

        Parameters
        ----------
        components : list
          List of AbsComponent objects
        """
        # Check
        assert ltiu.chk_components(components)
        # Instantiate with the first component
        init_comp = components[0]
        slf = cls(init_comp.coord, init_comp.zcomp, init_comp.vlim)
        if slf.chk_component(init_comp):
            slf._components.append(init_comp)
        else:
            raise IOError("Bad component input")
        # Append with component checking
        if len(components) > 1:
            for component in components[1:]:
                slf.add_component(component)
        # Return
        return slf

    def __init__(self, abs_type, radec, zabs, vlim, zem=0., NHI=0., sig_NHI=np.zeros(2), name=''):

        self.zabs = zabs
        self.zem = zem
        self.vlim = vlim
        self.NHI = NHI
        self.sig_NHI = sig_NHI
        # RA/DEC
        if isinstance(radec,(tuple)):
            self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        elif isinstance(radec,SkyCoord):
            self.coord = radec
        self.name = name

        # Abs type
        if abs_type == None:
            self.abs_type = 'NONE'
        else:
            self.abs_type = abs_type

        # Components
        self._components = []  # List of AbsComponent objects

        # Kinematics
        self.kin = {}

        # Abundances
        self._EW = QTable()
        self._ionN = QTable()
        self._ionstate = {}
        self._abund = QTable()

        # Refs (list of references)
        self.Refs = []

    def add_component(self, abscomp, toler=0.2*u.arcsec):
        """Add an AbsComponent object if it satisfies all of the rules.

        For velocities, we demand that the new component has a velocity
        range that is fully encompassed by the AbsSystem.

        Should check for duplicates..

        Parameters
        ----------
        comp : AbsComponent
        toler : Angle, optional
          Tolerance on matching coordinates
        """
        # Coordinates
        test = bool(self.coord.separation(abscomp.coord) < toler)
        # Now redshift/velocity
        zlim_comp = (1+abscomp.zcomp)*abscomp.vlim/const.c.to('km/s')
        zlim_sys = (1+self.zabs)*self.vlim/const.c.to('km/s')
        test = test & (zlim_comp[0]>=zlim_sys[0]) & (zlim_comp[1]<=zlim_sys[1])

        # Additional checks (specific to AbsSystem type)
        test = test & self.chk_component(abscomp)

        # Append?
        if test:
            self._components.append(abscomp)
        else:
            warnings.warn('Input AbsComponent with does not match AbsSystem rules. Not appending')

    # #############
    def __repr__(self):
        txt = '[{:s}: name={:s} type={:s}, {:s} {:s}, z={:g}, NHI={:g}'.format(
            self.__class__.__name__, self.name, self.abs_type,
            self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
            self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
            self.zabs, self.NHI)
        # Finish
        txt = txt + ']'
        return (txt)

    def chk_component(self, component):
        """Additional checks on the component"""
        return True

    def list_of_abslines(self):
        """ Generate a list of the absorption lines in this system

        Drawn from the components

        Returns
        -------
        abslist : list

        """
        abslist = []
        for component in self._components:
            for iline in component._abslines:
                abslist.append(iline)
        # Return
        return abslist


class GenericAbsSystem(AbsSystem):
    """Class for Generic Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
        if vlim is None:
            vlim = [-500.,500.]*u.km/u.s
        AbsSystem.__init__(self, 'Generic', radec, zabs, vlim, **kwargs)
        self.name = 'Foo'

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

class LymanAbsSystem(AbsSystem):
    """Class for HI Lyman Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
        AbsSystem.__init__(self, 'HILyman', radec, zabs, vlim, **kwargs)

    def chk_component(self,component):
        """Require components are only of HI
        """
        # Require HI
        test = (component.Zion[0] == 1) & (component.Zion[1] == 1)
        if not test:
            warnings.warn('Input component must be HI')
        return test

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'HILyman'

class AbsSubSystem(object):
    """
    Sub system.  Most frequently used in LLS

    Used to analyze a portion of an AbsSystem
    Must be tied to an AbsSystem

    Parameters
    ----------
    parent : AbsSystem
      Link
    zabs : float
      Absorption redshift
    vlim : Quantity array (2)
      Velocity limits of the system
    lbl : str
      Label for the SubSystem, e.g. 'A', 'B'
    """
    def __init__(self, parent, zabs, vlim, lbl):
        self.parent = parent
        self.zabs = zabs
        self.vlim = vlim
        self.lbl = lbl

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'SubSystem'

        # #############
    def __repr__(self):
        txt = '[{:s}: name={:s}{:s} type={:s}, {:s} {:s}, z={:g}, vlim={:g},{:g}'.format(
            self.__class__.__name__, self.parent.name, self.lbl,
            self.parent.abs_type,
            self.parent.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
            self.parent.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
            self.zabs, self.vlim[0],self.vlim[1])
        # Finish
        txt = txt + ']'
        return (txt)
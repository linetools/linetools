""" Class for  absorption systems
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import pdb
import warnings
from abc import ABCMeta

from astropy import units as u
from astropy.table import QTable
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu
from linetools import utils as ltu
from linetools.abund import ions


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

    @classmethod
    def from_dict(cls, idict):
        """ Instantiate from a dict.  Usually read from the hard-drive

        Parameters
        ----------
        idict : dict

        Returns
        -------
        AbsSystem

        """
        slf = cls(idict['abs_type'],
                  SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg),
                  idict['zabs'], idict['vlim']*u.km/u.s,
                  zem=idict['zem'], NHI=idict['NHI'],
                  sig_NHI=idict['sig_NHI'], name=idict['Name']
                  )
        # Components
        for key in idict['components']:
            comp = AbsComponent.from_dict(idict['components'][key])
            slf.add_component(comp)
        # Return
        return slf


    def __init__(self, abs_type, radec, zabs, vlim, zem=0.,
                 NHI=0., sig_NHI=np.zeros(2), name=None):

        self.zabs = zabs
        self.zem = zem
        self.vlim = vlim
        self.NHI = NHI
        self.sig_NHI = sig_NHI
        self.coord = ltu.radec_to_coord(radec)
        if name is None:
            self.name = 'J{:s}{:s}_z{:.3f}'.format(
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True),
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True),
                    self.zabs)
        else:
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

    def chk_component(self, component):
        """Additional checks on the component

        Parameters
        ----------
        component : AbsComponent
        """
        return True

    def list_of_abslines(self):
        """ Generate a list of the absorption lines in this system

        Drawn from the components

        Returns
        -------
        abslist : list of AbsLines

        """
        # Return
        return [iline for component in self._components
                for iline in component._abslines]

    def to_dict(self):
        """ Write AbsSystem data to a dict
        """
        outdict = dict(Name=self.name, abs_type=self.abs_type, zabs=self.zabs,
                       vlim=self.vlim.to('km/s').value, zem=self.zem,
                       NHI=self.NHI, sig_NHI=self.sig_NHI,
                       RA=self.coord.ra.value, DEC=self.coord.dec.value,
                       kin=self.kin, Refs=self.Refs,
                       )
        # Components
        outdict['components'] = {}
        for component in self._components:
            outdict['components'][component.name] = component.to_dict()
        # Polish
        outdict = ltu.jsonify_dict(outdict)
        # Return
        return outdict

    def __repr__(self):
        txt = '<{:s}: name={:s} type={:s}, {:s} {:s}, z={:g}, NHI={:g}'.format(
                self.__class__.__name__, self.name, self.abs_type,
                self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                self.zabs, self.NHI)
        # Finish
        txt = txt + '>'
        return (txt)


class GenericAbsSystem(AbsSystem):
    """Class for Generic Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
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


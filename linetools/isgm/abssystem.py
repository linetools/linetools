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
from astropy.units import Quantity
from astropy.table import QTable
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm import utils as ltiu
from linetools import utils as ltu
from linetools.spectralline import AbsLine
from linetools.abund import ions

# Globals to speed things up
c_mks = const.c.to('km/s')

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
    def from_abslines(cls, abslines, vlim=None, **kwargs):
        """Instantiate from a list of AbsLines

        Parameters
        ----------
        components : list
          List of AbsComponent objects
        """
        # Generate components
        components = ltiu.build_components_from_abslines(abslines, **kwargs)
        # Instantiate
        slf = cls.from_components(components, vlim=vlim)
        # Return
        return slf

    @classmethod
    def from_components(cls, components, vlim=None, NHI=None):
        """Instantiate from a list of AbsComponent objects

        Parameters
        ----------
        components : list
          List of AbsComponent objects
        vlim : list, optional
          Velocity limits for the system
          If not set, the first components sets vlim
        NHI : float, optional
          Set the NHI value of the system.  If not set,
          the method sums the NHI values of all the HI
          components input (if any)
        """
        # Check
        assert ltiu.chk_components(components)
        # Instantiate with the first component
        init_comp = components[0]
        if vlim is None:
            vlim = init_comp.vlim
        # Attempt to set NHI
        HI_comps = [comp for comp in components if comp.Zion == (1,1)]
        if NHI is None:
            NHI = 0.
            for HI_comp in HI_comps:  # Takes only the first line in each list
                NHI += HI_comp._abslines[0].attrib['N'].value
            if NHI > 0.:
                NHI = np.log10(NHI)
        #
        slf = cls(init_comp.coord, init_comp.zcomp, vlim, NHI=NHI)
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
    def from_json(cls, json_file, **kwargs):
        """ Load from a JSON file (via from_dict)

        Parameters
        ----------
        json_file
        kwargs

        Returns
        -------
        AbsSystem

        """
        idict = ltu.loadjson(json_file)
        slf = cls.from_dict(idict, **kwargs)
        return slf

    @classmethod
    def from_dict(cls, idict, skip_components=False, use_coord=False, **kwargs):
        """ Instantiate from a dict.  Usually read from the hard-drive

        Parameters
        ----------
        idict : dict
        skip_components : bool, optional
          If True, absorption components (if any exist) are not loaded from the input dict.
          Use when you are only interested in the global properties of an AbsSystem
        use_coord : bool, optinal
          Use coordinates from the AbsSystem to build the components (and lines)
          Speeds up performance, but you should know things are OK before using this

        Returns
        -------
        AbsSystem

        """
        #slf = cls(idict['abs_type'], SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg), idict['zabs'], idict['vlim']*u.km/u.s, zem=idict['zem'], NHI=idict['NHI'], sig_NHI=idict['sig_NHI'], flag_NHI=idict['flag_NHI'], name=idict['Name'] )
        slf = cls(SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg), idict['zabs'], idict['vlim']*u.km/u.s, zem=idict['zem'], NHI=idict['NHI'], sig_NHI=idict['sig_NHI'], flag_NHI=idict['flag_NHI'], name=idict['Name'] )
        if not skip_components:
            # Components
            if use_coord:  # Speed up performance
                coord = slf.coord
            else:
                coord = None
            components = ltiu.build_components_from_dict(idict, coord=coord, **kwargs)
            for component in components:
                # This is to insure the components follow the rules
                slf.add_component(component, **kwargs)

        # Return
        return slf

    def __init__(self, radec, zabs, vlim, zem=0., abs_type=None,
                 NHI=0., sig_NHI=np.zeros(2), flag_NHI=0, name=None):

        self.zabs = zabs
        self.zem = zem
        self.vlim = vlim
        self.NHI = NHI
        self.sig_NHI = sig_NHI
        self.flag_NHI = flag_NHI
        self.coord = ltu.radec_to_coord(radec)
        if name is None:
            self.name = 'J{:s}{:s}_z{:.3f}'.format(
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True),
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True),
                    self.zabs)
        else:
            self.name = name

        # Abs type
        if abs_type is None:
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

    def add_component(self, abscomp, tol=0.2*u.arcsec, chk_sep=True, **kwargs):
        """Add an AbsComponent object if it satisfies all of the rules.

        For velocities, we demand that the new component has a velocity
        range that is fully encompassed by the AbsSystem.

        Should check for duplicates..

        Parameters
        ----------
        comp : AbsComponent
        tol : Angle, optional
          Tolerance on matching coordinates
          Only used if chk_sep=True
        chk_sep : bool, optional
          Perform coordinate check (expensive)
        """
        # Coordinates
        if chk_sep:
            test = bool(self.coord.separation(abscomp.coord) < tol)
        else:
            test = True
        # Now redshift/velocity
        zlim_comp = (1+abscomp.zcomp)*abscomp.vlim/c_mks
        zlim_sys = (1+self.zabs)*self.vlim/c_mks
        test = test & (zlim_comp[0] >= zlim_sys[0]) & (zlim_comp[1] <= zlim_sys[1])

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

    def fill_ionN(self, **kwargs):
        """ Fills the ionN Table from the list of components
        """
        self._ionN = ltiu.iontable_from_components(self._components, **kwargs)

    def get_absline(self, inp):
        """ Returns an AbsLine from the AbsSystem

        Parameters
        ----------
        inp : str or Quantity
          str -- Name of the transition, e.g. 'CII 1334'
          Quantity -- Rest wavelength of the transition, e.g. 1334.53*u.AA
            to 0.01 precision

        Returns
        -------
        absline -- AbsLine object or list of Abslines
          More than one will be returned if this line exists in
          multiple components.  The returned quantity will then
          be a list instead of a single object
        """
        # Generate the lines
        abslines = self.list_of_abslines()
        if isinstance(inp,basestring):
            names = np.array([absline.name for absline in abslines])
            mt = np.where(names == inp)[0]
        elif isinstance(inp,Quantity):
            wrest = Quantity([absline.wrest for absline in abslines])
            mt = np.where(np.abs(wrest-inp) < 0.01*u.AA)[0]
        else:
            raise IOError("Bad input to absline")
        # Finish
        if len(mt) == 0:
            warnings.warn("No absline with input={}".format(inp))
            return None
        elif len(mt) == 1:
            return abslines[mt]
        else:
            return [abslines[ii] for ii in mt]

    def get_comp_from_absline(self, aline):
        """ Returns the component that holds the input AbsLine

        Parameters
        ----------
        aline : AbsLine

        Returns
        -------
        comp -- AbsComponent object that holds this AbsLine
        """
        # Loop on components
        for comp in self._components:
            # Is the line present?
            try:
                idx = comp._abslines.index(aline)
            except ValueError:
                pass
            else:
                return comp
        # Raise error?
        warnings.warn("Input absorption line is not in any component")
        return None

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
        """ Write AbsSystem data to a dict that can be written with JSON
        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        outdict = dict(Name=self.name, abs_type=self.abs_type, zabs=self.zabs,
                       vlim=self.vlim.to('km/s').value, zem=self.zem,
                       NHI=self.NHI, sig_NHI=self.sig_NHI, flag_NHI=self.flag_NHI,
                       RA=self.coord.ra.value, DEC=self.coord.dec.value,
                       kin=self.kin, Refs=self.Refs, CreationDate=date,
                       user=user
                       )
        # Components
        outdict['components'] = {}
        for component in self._components:
            outdict['components'][component.name] = ltu.jsonify(component.to_dict())
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

    def write_json(self, outfil=None):
        """ Generate a JSON file from the system

        Returns
        -------

        """
        import io, json
        # Generate the dict
        odict = self.to_dict()
        # Write
        if outfil is None:
            outfil = self.name+'.json'
        with io.open(outfil, 'w', encoding='utf-8') as f:
            f.write(json.dumps(odict, sort_keys=True, indent=4,
                               separators=(',', ': ')))
        # Finish
        print("Wrote {:s} system to {:s} file".format(self.name, outfil))


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
        AbsSystem.__init__(self, radec, zabs, vlim, abs_type='Generic', **kwargs)
        self.name = 'Foo'

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

class LymanAbsSystem(AbsSystem):
    """Class for HI Lyman Absorption Line System
    """
    def __init__(self, radec, zabs, vlim, **kwargs):
        AbsSystem.__init__(self, radec, zabs, vlim, abs_type='HILyman', **kwargs)

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


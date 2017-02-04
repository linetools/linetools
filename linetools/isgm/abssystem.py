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
from linetools import line_utils as ltlu
from linetools.spectralline import AbsLine
from linetools.abund import ions

# Globals to speed things up
c_mks = const.c.to('km/s').value

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
        ZH :  float
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
            for HI_comp in HI_comps:
                NHI += 10**HI_comp.logN
                #NHI += HI_comp._abslines[0].attrib['N'].value
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
    def from_dict(cls, idict, use_coord=False, **kwargs):
        """ Instantiate from a dict.  Usually read from the hard-drive

        Parameters
        ----------
        idict : dict

        Returns
        -------
        AbsSystem

        """
        if 'NHI' in idict.keys():
            ckwargs = dict(NHI=idict['NHI'], sig_NHI=idict['sig_NHI'], flag_NHI=idict['flag_NHI'])
        slf = cls(SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg),
                  idict['zabs'], idict['vlim']*u.km/u.s, zem=idict['zem'],
                  name=idict['Name'], **ckwargs)
        # Other
        add_other_from_dict(slf, idict)
        # Components
        add_comps_from_dict(slf, idict, **kwargs)

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

        # Metallicity
        self.ZH = 0.
        self.sig_ZH = 0.

        # Abundances and Tables
        self._EW = QTable()
        self._ionN = None   # Needs to be None for fill_ion
        self._trans = QTable()
        self._ionstate = {}
        self._abund = QTable()

        # Refs (list of references)
        self.Refs = []

    def add_component(self, abscomp, tol=0.2*u.arcsec,
                      chk_sep=True, chk_z=True, overlap_only=False,
                      vtoler=1., debug=False, **kwargs):
        """Add an AbsComponent object if it satisfies all of the rules.

        For velocities, we demand that the new component has a velocity
        range that is fully encompassed by the AbsSystem.
        We allow a small tolerance for round-off error

        Should check for duplicates..

        Parameters
        ----------
        comp : AbsComponent
        tol : Angle, optional
          Tolerance on matching coordinates
          Only used if chk_sep=True
        chk_sep : bool, optional
          Perform coordinate check (expensive)
        chk_z : bool, optional
          Perform standard velocity range test
        overlap_only : bool, optional
          Only require that the components overlap in redshift
        vtoler : float, optional
          Tolerance for velocity in km/s

        Returns
        -------
        test : bool
          True if successful
        """
        # Coordinates
        if chk_sep:
            testcoord = bool(self.coord.separation(abscomp.coord) < tol)
        else:
            testcoord = True
        # Now redshift/velocity
        testz = True
        if chk_z:
            # Will avoid Quantity for speed
            comp_vlim_mks = abscomp.vlim.to('km/s').value
            sys_vlim_mks = self.vlim.to('km/s').value
            dz_toler = (1 + self.zabs) * vtoler / c_mks
            zlim_comp = abscomp.zcomp + (1 + abscomp.zcomp) * (comp_vlim_mks / c_mks)
            zlim_sys = self.zabs + (1 + self.zabs) * (sys_vlim_mks / c_mks)
            if overlap_only:
                testz = True
                if debug:
                    pdb.set_trace()
                if np.all(zlim_comp > np.max(zlim_sys + dz_toler)) or np.all(
                                zlim_comp < np.min(zlim_sys-dz_toler)):
                    testz = False
            else:
                testz = (zlim_comp[0] >= (zlim_sys[0]-dz_toler)) & (
                    zlim_comp[1] <= (zlim_sys[1]+dz_toler))

        # Additional checks (specific to AbsSystem type)
        testcomp = self.chk_component(abscomp)
        test = testcoord & testcomp & testz

        # Append?
        if test:
            self._components.append(abscomp)
        else:
            warnings.warn('Input AbsComponent with Zion={} does not match AbsSystem rules. Not appending'.format(abscomp.Zion))
            if not testcoord:
                warnings.warn('Failed coordinate match')
            if not testcomp:
                warnings.warn('Failed component check')
            if not testz:
                warnings.warn('Failed velocity overlap')
        #
        return test

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

    def fill_trans(self, **kwargs):
        """ Fills the ionN Table from the list of components
        """
        self._trans = ltlu.transtable_from_speclines(self.list_of_abslines())

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
            return abslines[mt[0]]
        else:
            return [abslines[ii] for ii in mt]

    def get_component(self, inp):
        """ Returns the component related to the given input
        TODO: Need to handle fine-structure lines at some point..

        Parameters
        ----------
        inp : tuple or AbsLine
          tuple -- (Z,ion) integers
          AbsLine -- actual absorption line object

        Returns
        -------
        component
        """
        if isinstance(inp, tuple):
            # Assume Zion for now, e.g. (26,2)
            tuples = [comp.Zion for comp in self._components]
            try:
                idx = tuples.index(inp)
            except ValueError:
                warnings.warn("Input Zion {} is not in any component".format(inp))
                return None
            else:
                return self._components[idx]
        elif isinstance(inp, AbsLine):
            return self.get_comp_from_absline(inp)
        else:
            raise IOError("Bad input to get_component method")

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

    def measure_restew(self, spec=None, **kwargs):
        """ Measure rest-frame EWs for lines in the AbsSystem
        Analysis is only performed on lines with analy['do_analysis'] != 0

        Parameters
        ----------
        spec : XSpectrum1D, optional
        kwargs

        Returns
        -------

        """
        # Grab Lines
        abs_lines = self.list_of_abslines()
        # Loop
        for iline in abs_lines:
            # Fill in spec?
            if spec is not None:
                iline.analy['spec'] = spec
            # Check for analysis
            if iline.analy['do_analysis'] == 0:
                warnings.warn("Skipping {:s} because do_analysis=0".format(iline.name))
            else:
                # Measure
                iline.measure_restew(**kwargs)

    def measure_aodm(self, spec=None, **kwargs):
        """ Measure ADOM columns for the list of lines
        Note: Components are *not* updated by default

        Parameters
        ----------
        spec : XSpectrum1D, optional
        kwargs

        Returns
        -------

        """
        # Grab Lines
        abs_lines = self.list_of_abslines()
        # Loop
        for iline in abs_lines:
            # Fill in spec?
            if spec is not None:
                iline.analy['spec'] = spec
            # Measure
            iline.measure_aodm(**kwargs)
        #
        print("You may now wish to update the component column densities with update_component_colm()")

    def update_component_colm(self, **kwargs):
        """ Synthesize/update column density measurements for components

        Parameters
        ----------
        kwargs

        Returns
        -------

        """
        for comp in self._components:
            comp.synthesize_colm(**kwargs)

    def stack_plot(self, pvlim=None, **kwargs):
        """Show a stack plot of the system, if spec are loaded
        Assumes the data are normalized.

        Parameters
        ----------
        pvlim : Quantities, optional
          Over-ride system vlim for plotting
        """
        from linetools.analysis import plots as ltap
        if pvlim is not None:
            vlim = pvlim
        else:
            vlim = self.vlim
        ltap.stack_plot(self.list_of_abslines(), vlim=vlim, **kwargs)

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
                       ZH=self.ZH, sig_ZH=self.sig_ZH,
                       user=user
                       )
        outdict['class'] = self.__class__.__name__
        # Components
        outdict['components'] = {}
        for component in self._components:
            outdict['components'][component.name] = ltu.jsonify(component.to_dict())
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

    def update_vlim(self, sub_system=None):
        """ Update vlim in the main or subsystems

        Parameters
        ----------
        sub_system : str, optional
          If provided, apply to given sub-system.  Only used in LLS so far
        """
        def get_vmnx(components):
            zlim_sys = ltu.z_from_dv(self.vlim, self.zabs, rel=False)
            zmin, zmax = zlim_sys
            for component in components:
                zlim_comp = ltu.z_from_dv(component.vlim, component.zcomp, rel=False)
                zmin = min(zmin, zlim_comp[0])
                zmax = max(zmax, zlim_comp[1])
            # Convert back to velocities
            return ltu.dv_from_z([zmin,zmax], self.zabs, rel=False)

        # Sub-system?
        if sub_system is not None:
            components = self.subsys[sub_system]._components
            self.subsys[sub_system].vlim = get_vmnx(components)
        else:
            components = self._components
            self.vlim = get_vmnx(components)  # Using system z

    def write_json(self, outfil=None):
        """ Generate a JSON file from the system

        Returns
        -------

        """
        # Generate the dict
        odict = self.to_dict()
        # Write
        if outfil is None:
            outfil = self.name+'.json'
        ltu.savejson(outfil, odict, overwrite=True, easy_to_read=True)
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


def add_comps_from_dict(slf, idict, skip_components=False, use_coord=False, **kwargs):
    """
    Parameters
    ----------
    slf : AbsSystem, AbsSightline
      Or any object with an add_component() method
    skip_components : bool, optional
      If True, absorption components (if any exist) are not loaded from the input dict.
      Use when you are only interested in the global properties of an AbsSystem
    use_coord : bool, optinal
      Use coordinates from the AbsSystem to build the components (and lines)
      Speeds up performance, but you should know things are OK before using this

    Returns
    -------

    """
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

def add_other_from_dict(slf, idict):
    """ Add other attributes to a system from a dict
    Useful for handling various AbsSystem types

    Parameters
    ----------
    slf : AbsSystem
    idict : dict
    """
    # Other
    if 'kin' in idict.keys():
        slf.kin = ltu.convert_quantity_in_dict(idict['kin'])
    if 'Refs' in idict.keys():
        slf.Refs = idict['Refs']
    if 'ZH' in idict.keys():
        slf.ZH = idict['ZH']
        slf.sig_ZH = idict['sig_ZH']

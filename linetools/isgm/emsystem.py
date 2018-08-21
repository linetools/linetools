""" Class for  emission line systems
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
from astropy.table import Table
from astropy import constants as const
from astropy.coordinates import SkyCoord

from linetools.isgm import utils as ltiu
from linetools import utils as ltu
from linetools import line_utils as ltlu
from linetools.spectralline import EmLine
from linetools import io as lio

# Globals to speed things up
c_mks = const.c.to('km/s').value

class EmSystem(object):
    """
    Class for an emission line system, e.g. Planetary Nebula, Galaxy, SN, quasar

    Parameters
    ----------
    radec : tuple or coordinate
        RA/Dec of the sightline or astropy.coordinate
    zem : float
      Emission redshift
    vlim : Quantity array (2)
      Velocity limits of the system
    zem : float, optional
      Emission redshift of background source
    em_type : str or unicode
      Type of emission system, e.g. galaxy, quasar, LAE
    name : str, optional
      Name for the system

    Attributes
    ----------
        em_type : str
        coord : SkyCoord
            RA/Dec of the sightline
        zem : float
          Emission redshift
        zem : float
          Emission redshift of background source
        vlim : Quantity array (2)
          Velocity limits of the system
        ZH :  float
          Metallicity (log10)
        name : str
            Name of the system
    """

    __metaclass__ = ABCMeta

    @classmethod
    def from_emlines(cls, emlines, vlim=None, **kwargs):
        """Instantiate from a list of EmLines

        Parameters
        ----------
        emlines : list
          List of EmLine objects
        vlim : list, optional
          Velocity limits for the system
          If not set, the first components sets vlim
        """
        # Check
        #assert ltiu.chk_components(components)
        # Instantiate with the first component
        init_eml = emlines[0]
        #
        slf = cls(init_eml.attrib['coord'], init_eml.z, vlim=vlim)
        if slf.chk_emline(init_eml):
            slf._emlines.append(init_eml)
        else:
            raise IOError("Bad component input")
        # Append with component checking
        for emline in emlines[1:]:
            slf.add_emline(emline)
        # Return
        return slf

    @classmethod
    def from_alis(cls, alis_file, radec, **kwargs):
        """ Load from an ALIS output file

        Parameters
        ----------
        alis_file : .mod ALIS output file
        radec : SkyCoord or tuple of RA,DEC
        kwargs

        Returns
        -------
        EmSystem

        """
        emlines = lio.emlines_from_alis_output(alis_file)
        # Add coordinates
        coord = ltu.radec_to_coord(radec)
        for emline in emlines:
            emline.attrib['coord'] = coord
        #
        return cls.from_emlines(emlines, **kwargs)

    @classmethod
    def from_dict(cls, idict, use_coord=False, **kwargs):
        """ Instantiate from a dict.  Usually read from the hard-drive

        Parameters
        ----------
        idict : dict

        Returns
        -------

        """
        slf = cls(SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg),
                  idict['zem'], vlim=idict['vlim']*u.km/u.s,
                  name=idict['Name'], em_type=idict['em_type'])
        # Emission lines
        em_dict = idict['emlines']
        for key in em_dict:
            iline = em_dict[key]
            obj = EmLine.from_dict(iline)
            slf.add_emline(obj)
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
        EmSystem

        """
        idict = ltu.loadjson(json_file)
        slf = cls.from_dict(idict, **kwargs)
        return slf


    def __init__(self, radec, zem, vlim=None, em_type=None, name=None):

        self.zem = zem
        if vlim is None:
            self.vlim = [-300., 300.]*u.km/u.s
        else:
            self.vlim = vlim
        self.coord = ltu.radec_to_coord(radec)
        if name is None:
            self.name = 'J{:s}{:s}_z{:.3f}'.format(
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True),
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True),
                    self.zem)
        else:
            self.name = name

        # Em type
        if em_type is None:
            self.em_type = 'NONE'
        else:
            self.em_type = em_type

        # Components
        self._emlines = []  # List of EmLine objects

        # Components
        #self._components = []  # List of EmComponent objects

        # Kinematics
        self.kin = {}

        # Metallicity
        self.ZH = 0.
        self.sig_ZH = 0.

        # Abundances and Tables
        self._EW = Table()
        self._fluxes = None   # Needs to be None for fill_ion
        #self._trans = Table()
        self._abund = Table()

        # Refs (list of references)
        self.Refs = []

    def add_emlines_from_alis(self, alis_file, **kwargs):
        '''
        Code assumes the coordinates of the class (i.e. no checking)
        Parameters
        ----------
        alis_file : str



        Returns
        -------

        '''
        dum_sys = EmSystem.from_alis(alis_file, self.coord)
        # Add in lines
        for emline in dum_sys._emlines:
            self.add_emline(emline, **kwargs)

    def add_emline(self, emline, tol=0.2*u.arcsec,
                      chk_sep=True, chk_z=True, overlap_only=False,
                      vtoler=1., debug=False, **kwargs):
        """Add an EmLine object if it satisfies all of the rules.

        For velocities, we demand that the new line has a velocity
        range that is fully encompassed by the EmSystem.
        We allow a small tolerance for round-off error

        Should check for duplicates..

        Parameters
        ----------
        emline : EmLine
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
            testcoord = bool(self.coord.separation(emline.attrib['coord']) < tol)
        else:
            testcoord = True
        # Now redshift/velocity
        if chk_z:
            # Will avoid Quantity for speed
            sys_vlim_mks = self.vlim.to('km/s').value
            dz_toler = (1 + self.zem) * vtoler / c_mks
            zlim_sys = self.zem + (1 + self.zem) * (sys_vlim_mks / c_mks)
            testz = (emline.z >= (zlim_sys[0]-dz_toler)) & (
                    emline.z <= (zlim_sys[1]+dz_toler))
        else:
            testz = True

        # Additional checks (specific to EmSystem type)
        testcomp = self.chk_emline(emline)
        test = testcoord & testcomp & testz

        # Append?
        if test:
            self._emlines.append(emline)
        else:
            warnings.warn('Input EmLine with wrest={} does not match EmSystem rules. Not appending'.format(emline.wrest))
            if not testcoord:
                warnings.warn('Failed coordinate match')
            if not testcomp:
                warnings.warn('Failed component check')
            if not testz:
                warnings.warn('Failed velocity overlap')
        #
        return test

    def chk_emline(self, component):
        """Additional checks on the component

        Parameters
        ----------
        component : EmComponent

        """
        return True

    def fill_ionN(self, **kwargs):
        """ Fills the ionN Table from the list of components
        """
        self._ionN = ltiu.iontable_from_components(self._components, **kwargs)

#    def fill_trans(self, **kwargs):
#        """ Fills the ionN Table from the list of emission lines
#        """
#        self._trans = ltlu.transtable_from_speclines(self.list_of_abslines())

    def get_emline(self, inp):
        """ Returns an EmLine from the EmSystem

        Parameters
        ----------
        inp : str or Quantity
          str -- Name of the transition, e.g. 'Halpha'
          Quantity -- Rest wavelength of the transition, e.g. 6564.61*u.AA
            to 0.01 precision

        Returns
        -------
        emline -- EmLine object or list of EmLines
          More than one will be returned if this line exists in
          multiple components.  The returned quantity will then
          be a list instead of a single object
        """
        # Generate the lines
        if isinstance(inp,basestring):
            names = np.array([emline.name for emline in self._emlines])
            mt = np.where(names == inp)[0]
        elif isinstance(inp,Quantity):
            wrest = Quantity([emline.wrest for emline in self._emlines])
            mt = np.where(np.abs(wrest-inp) < 0.01*u.AA)[0]
        else:
            raise IOError("Bad input to emline")
        # Finish
        if len(mt) == 0:
            warnings.warn("No emline with input={}".format(inp))
            return None
        elif len(mt) == 1:
            return self._emlines[mt[0]]
        else:
            return [self._emlines[ii] for ii in mt]

#    def get_component(self, inp):
#        """ Returns the component related to the given input
#        Parameters
#        ----------
#        inp : tuple or AbsLine
#          tuple -- (Z,ion) integers
#          AbsLine -- actual absorption line object

#        Returns
#        -------
#        component
#        """
#        if isinstance(inp, tuple):
#            # Assume Zion for now, e.g. (26,2)
#            tuples = [comp.Zion for comp in self._components]
#            try:
#                idx = tuples.index(inp)
#            except ValueError:
#                warnings.warn("Input Zion {} is not in any component".format(inp))
#                return None
#            else:
#                return self._components[idx]
#        elif isinstance(inp, AbsLine):
#            return self.get_comp_from_absline(inp)
#        else:
#            raise IOError("Bad input to get_component method")

#    def get_comp_from_absline(self, aline):
#        """ Returns the component that holds the input AbsLine

#        Parameters
#        ----------
#        aline : AbsLine

#        Returns
#        -------
#        comp -- AbsComponent object that holds this AbsLine
#        """
#        # Loop on components
#        for comp in self._components:
#            # Is the line present?
#            try:
#                idx = comp._abslines.index(aline)
#            except ValueError:
#                pass
#            else:
#                return comp
#        # Raise error?
#        warnings.warn("Input absorption line is not in any component")
#        return None

    def measure_restew(self, spec=None, **kwargs):
        """ Measure rest-frame EWs for lines in the EmSystem
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
        fig = ltap.stack_plot(self.list_of_abslines(), vlim=vlim, **kwargs)
        if fig is not None:
            return fig

    def to_dict(self):
        """ Write EmSystem data to a dict that can be written with JSON
        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        outdict = dict(Name=self.name, em_type=self.em_type,
                       vlim=self.vlim.to('km/s').value, zem=self.zem,
                       RA=self.coord.ra.value, DEC=self.coord.dec.value,
                       kin=self.kin, Refs=self.Refs, CreationDate=date,
                       ZH=self.ZH, sig_ZH=self.sig_ZH,
                       user=user
                       )
        outdict['class'] = self.__class__.__name__
        outdict['emlines'] = {}
        for iline in self._emlines:
            outdict['emlines'][iline.wrest.value] = iline.to_dict()
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

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
        txt = '<{:s}: name={:s} type={:s}, {:s} {:s}, z={:g}'.format(
                self.__class__.__name__, self.name, self.em_type,
                self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                self.zem)
        # Finish
        txt = txt + '>'
        return (txt)


class GenericEmSystem(EmSystem):
    """Class for Generic Emission Line System
    """
    def __init__(self, radec, zem, **kwargs):
        EmSystem.__init__(self, radec, zem, em_type='Generic', **kwargs)
        self.name = 'Foo'

    def print_em_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'Generic'

""" Classes for absorption line component
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb
import numpy as np
import warnings
from abc import ABCMeta

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column

from linetools.analysis import absline as ltaa
from linetools.spectralline import AbsLine, SpectralLine
from linetools import utils as ltu

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# Global import for speed
c_kms = const.c.to('km/s').value

# Class for Sightline
class AbsSightline(object):
    """ Abstract Class for an absorption sightline

    Attributes
    ----------
    name : str
        Name of the component, e.g. `Si II`
    coord : SkyCoord
        Sky coordinate
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
    comment : str
        A comment, default is ``
    """
    __metaclass__ = ABCMeta

    @classmethod
    def from_abslines(cls, abslines, stars=None, **kwargs):
        """Instantiate from a list of AbsLine objects

        Parameters
        ----------
        abslines : list 
          List of AbsLine objects
        stars : str, optional
          Asterisks to append to the ion name (e.g. fine-structure, CII*)
        """
        # Check
        if not isinstance(abslines, list):
            raise IOError("Need a list of AbsLine objects")
        if not all(isinstance(x, AbsLine) for x in abslines):
            raise IOError("List needs to contain only AbsLine objects")

        # Instantiate with the first line
        init_line = abslines[0]
        #init_line.attrib['z'], init_line.analy['vlim'],
        slf = cls( init_line.attrib['coord'], (init_line.data['Z'],init_line.data['ion']),
                   init_line.attrib['z'], init_line.limits.vlim,
                   Ej=init_line.data['Ej'], stars=stars)
        slf._abslines.append(init_line)
        # Append with component checking
        if len(abslines) > 1:
            for absline in abslines[1:]:
                slf.add_absline(absline, **kwargs)
        # Return
        return slf

    @classmethod
    def from_component(cls, component, **kwargs):
        """ Instantiate from an AbsComponent object

        Uses RA/DEC, Zion, Ej, A, z, vlim

        Parameters
        ----------
        component : AbsComponent
           An AbsComponent object

        Returns
        -------
        AbsComponent
        """
        # Check
        if not isinstance(component, AbsComponent):
            raise IOError('Need an AbsComponent object')
        # Return
        return cls(component.coord, component.Zion, component.zcomp, component.vlim, Ej=component.Ej,
                   A=component.A, name=component.name, **kwargs)

    @classmethod
    def from_dict(cls, idict, coord=None, **kwargs):
        """ Instantiate from a dict

        Parameters
        ----------
        idict : dict

        Returns
        -------

        """
        if coord is not None:
            radec = coord
        else:
            radec = SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg)
        # Init
        #slf = cls(radec, tuple(idict['Zion']), idict['zcomp'], Quantity(idict['vlim'], unit='km/s'),
        slf = cls(radec, tuple(idict['Zion']), idict['zcomp'], idict['vlim']*u.km/u.s,
                  Ej=idict['Ej']/u.cm, A=idict['A'],
                  Ntup = tuple([idict[key] for key in ['flag_N', 'logN', 'sig_logN']]),
                  comment=idict['comment'], name=idict['Name'])
        # Add lines
        for key in idict['lines'].keys():
            iline = SpectralLine.from_dict(idict['lines'][key], coord=coord, **kwargs)
            slf.add_absline(iline, **kwargs)
        # Return
        return slf

    def __init__(self, radec, sl_type=None, em_type=None, comment=None, name=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
            (RA,DEC) in deg or astropy.coordinate
        sl_type : str, optional
          Sightline type, e.g. IGM
        emtype : str, optional
          Type of emission source for absorption sightline
        comment : str, optional
          A comment, default is ``
        """
        # Required
        self.coord = ltu.radec_to_coord(radec)

        # Other
        self._components = []

        # Name
        if name is None:
            self.name = 'J{:s}{:s}'.format(
                    self.coord.ra.to_string(unit=u.hour,sep='',pad=True),
                    self.coord.dec.to_string(sep='',pad=True,alwayssign=True))
        else:
            self.name=name

        # Others
        self.em_type = em_type
        self.sl_type = sl_type

    def add_absline(self, absline, tol=0.1*u.arcsec, chk_vel=True,
                    chk_sep=True, vtoler=1., **kwargs):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.

        For velocities, we demand that the new line has a velocity
        range that is fully encompassed by the component.

        Parameters
        ----------
        absline : AbsLine
        tol : Angle, optional
          Tolerance on matching coordinates.  Only used if chk_sep=True
        chk_vel : bool, optional
          Perform velocity test (can often be skipped)
          Insist the bounds of the AbsLine are within 1km/s of the Component
             (allows for round-off error)
        chk_sep : bool, optional
          Perform coordinate check (expensive)
        vtoler : float
          Tolerance for velocity in km/s (must be positive)
        """
        if vtoler < 0:
            raise ValueError('vtoler must be positive!')

        # Perform easy checks
        if chk_sep:
            testc = bool(self.coord.separation(absline.attrib['coord']) < tol)
        else:
            testc = True
        testZ = self.Zion[0] == absline.data['Z']
        testi = self.Zion[1] == absline.data['ion']
        testE = bool(self.Ej == absline.data['Ej'])
        # Now redshift/velocity
        if chk_vel:
            dz_toler = (1 + self.zcomp) * vtoler / c_kms  # Avoid Quantity for speed
            zlim_line = (1 + absline.attrib['z']) * absline.limits.vlim.to('km/s').value / c_kms
            zlim_comp = (1+self.zcomp) * self.vlim.to('km/s').value / c_kms
            testv = (zlim_line[0] >= (zlim_comp[0] - dz_toler)) & (
                zlim_line[1] <= (zlim_comp[1] + dz_toler))
        else:
            testv = True
        # Combine
        test = testc & testZ & testi & testE & testv
        # Isotope
        if self.A is not None:
            raise ValueError('Not ready for this yet.')
        # Append?
        if test:
            self._abslines.append(absline)
        else:
            warnings.warn("Failed add_absline test")
            print('Input absline with wrest={:g} does not match component rules. Not appending'.format(absline.wrest))
            if not testv:
                print("Absline velocities lie beyond component\n Set chk_vel=False to skip this test.")
            if not testc:
                print("Absline coordinates do not match.  Best to set them")

    def build_table(self):
        """Generate an astropy QTable out of the component.
        Returns
        -------
        comp_tbl : QTable
        """
        if len(self._abslines) == 0:
            return
        comp_tbl = QTable()
        comp_tbl.add_column(Column([iline.wrest.to(u.AA).value for iline in self._abslines]*u.AA, name='wrest'))
        for attrib in ['z', 'flag_N', 'logN', 'sig_logN']:
            comp_tbl.add_column(Column([iline.attrib[attrib] for iline in self._abslines], name=attrib))
        # Return
        return comp_tbl

    def to_dict(self):
        """ Write AbsSightline data to a dict that can be written with JSON
        """
        import datetime
        import getpass
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        user = getpass.getuser()
        # Generate the dict
        outdict = dict(RA=self.coord.ra.value, DEC=self.coord.dec.value,
                       CreationDate=date, user=user)
        # Add other attributes
        all_attr = self.__dict__
        for key in all_attr:
            if key in ['coord', '_components']:
                continue
            else:
                outdict[key] = getattr(self, key)
        # Components
        outdict['components'] = {}
        for component in self._components:
            outdict['components'][component.name] = ltu.jsonify(component.to_dict())
        # Polish
        outdict = ltu.jsonify(outdict)
        # Return
        return outdict

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}'.format(
            self.__class__.__name__, self.coord.ra.to_string(unit=u.hour,sep=':', pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True))

        # Type?
        if self.em_type is not None:
            txt = txt + ', emtype={:s}'.format(self.em_type)

        # Finish
        txt += '>'
        return (txt)


class GenericSightline(AbsSightline):
    """Class for Generic Absorption Sightline
    """
    def __init__(self, radec, **kwargs):
        AbsSightline.__init__(self, radec, sl_type='Generic', **kwargs)


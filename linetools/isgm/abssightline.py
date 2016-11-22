""" Class for an absorption sightline
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
from astropy.table import Table, Column

from linetools.spectralline import AbsLine, SpectralLine
from linetools import utils as ltu
from linetools.isgm.abscomponent import AbsComponent

# Global import for speed
c_kms = const.c.to('km/s').value

# Class for Sightline
class AbsSightline(object):
    """ Abstract Class for an absorption sightline

    Attributes
    ----------
    radec : SkyCoord or ra/dec tuple
        Coordinates
    sl_type : str, optional
        Description of type of sightline (e.g. IGM)
    em_type : str, optional
        Description of type of emission source (e.g. QSO, GRB)
    comment : str, optional
        A comment, default is ``
    name : str, optional
        Name of the sightline, e.g. '3C273'
    """
    __metaclass__ = ABCMeta

    @classmethod
    def from_abslines(cls, abslines, **kwargs):
        """Instantiate from a list of AbsLine objects

        Parameters
        ----------
        abslines : list 
          List of AbsLine objects
        stars : str, optional
          Asterisks to append to the ion name (e.g. fine-structure, CII*)
        """
        from .utils import build_components_from_abslines
        # Check
        if not isinstance(abslines, list):
            raise IOError("Need a list of AbsLine objects")
        if not all(isinstance(x, AbsLine) for x in abslines):
            raise IOError("List needs to contain only AbsLine objects")

        # Generate components
        comps = build_components_from_abslines(abslines, **kwargs)
        # Generate the sightline
        slf = cls.from_components(comps, **kwargs)

        # Return
        return slf

    @classmethod
    def from_components(cls, components, **kwargs):
        """ Instantiate from a list of AbsComponent objects

        Uses RA/DEC

        Parameters
        ----------
        components : list
           list of AbsComponent objects

        Returns
        -------
        """
        # Check
        if not isinstance(components, list):
            raise IOError("Need a list of AbsComponent objects")
        if not isinstance(components[0], AbsComponent):
            raise IOError('Need an AbsComponent object')

        # Instantiate with the first component
        slf = cls( components[0].coord)
        for comp in components:
            slf.add_component(comp, **kwargs)
        # Return
        return slf

    def __init__(self, radec, sl_type=None, em_type=None, comment=None, name=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
            (RA,DEC) in deg or astropy.coordinate.SkyCoord
        sl_type : str, optional
          Sightline type, e.g. IGM
        em_type : str, optional
          Type of emission source for absorption sightline, e.g. QSO
        comment : str, optional
          A comment, default is ``
        name : str, optional
            Name of the sightline, e.g. '3C273'
        """
        # Required
        self.coord = ltu.radec_to_coord(radec)

        # Lists
        self._components = []

        # Name
        if name is None:
            self.name = ltu.name_from_coord(self.coord)
        else:
            self.name = name

        # Others
        self.em_type = em_type
        self.sl_type = sl_type
        self._abssystems = None  # Saving the namespace for future usage

    def add_component(self, abscomp, tol=0.2*u.arcsec,
                      chk_sep=True, debug=False, **kwargs):
        """Add a component to AbsSightline if it satisfies all of the rules.

        Presently, the only constraint is on RA/DEC
        Note: It is likely that this method will be over-rode by the
        sub-class, e.g. to include source redshift in the criteria.

        Parameters
        ----------
        abscomp : AbsComponent
        tol : Angle, optional
          Tolerance on matching coordinates
          Only used if chk_sep=True
        chk_sep : bool, optional
          Perform coordinate check (expensive)
        """
        # Coordinates
        if chk_sep:
            testcoord = bool(self.coord.separation(abscomp.coord) < tol)
        else:
            testcoord = True

        # Combine (with other tests, if they exist)
        test = testcoord
        # Append?
        if test:
            self._components.append(abscomp)
        else:
            warnings.warn("Failed add_component test")
            #print('Input absline with wrest={:g} does not match component rules. Not appending'.format(absline.wrest))
            if not testcoord:
                print("AbsComp coordinates do not match.  Best to set them")

    def build_table(self):
        """Generate an astropy QTable out of the components.
        Default columns are z, Ion, flag_N, logN, sig_logN

        Returns
        -------
        comp_tbl : Table
        """
        from linetools.isgm.utils import iontable_from_components
        comp_tbl = Table()
        if len(self._components) == 0:
            return comp_tbl
        comp_tbl = iontable_from_components(self._components)

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
        outdict['class'] = self.__class__.__name__
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

        # Name
        if self.name is not None:
            txt = txt + ', name={:s}'.format(self.name)
        # Type?
        if self.em_type is not None:
            txt = txt + ', emtype={:s}'.format(self.em_type)

        # Finish
        txt += '>'
        return (txt)


class GenericAbsSightline(AbsSightline):
    """Class for Generic Absorption Sightline
    """
    def __init__(self, radec, **kwargs):
        AbsSightline.__init__(self, radec, sl_type='Generic', **kwargs)


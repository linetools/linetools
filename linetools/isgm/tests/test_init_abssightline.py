# Module to run tests on generating AbsSightline

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import GenericAbsSightline
from linetools.isgm.abssystem import GenericAbsSystem
from .utils import lyman_comp, si2_comp


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_from_systems():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    lyacomp1 = lyman_comp(radec, z=0.45)
    lyacomp2 = lyman_comp(radec, z=0.25)
    # SiII
    SiII_comp1 = si2_comp(radec, z=0.45)
    SiII_comp2 = si2_comp(radec, z=0.25)
    # Instantiate systems
    sys1 = GenericAbsSystem.from_components([lyacomp1,SiII_comp1])
    sys2 = GenericAbsSystem.from_components([lyacomp2, SiII_comp2])
    # Instantiate sightline
    sl = GenericAbsSightline.from_systems([sys1,sys2])
    # Test
    assert len(sl._abssystems) == 2

def test_from_components():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    # Init
    gensl = GenericAbsSightline.from_components([abscomp, SiII_comp])
    assert len(gensl._components) == 2


def test_from_abslines():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    # Init
    abs_lines = abscomp._abslines + SiII_comp._abslines
    gensl = GenericAbsSightline.from_abslines(abs_lines)
    assert len(gensl._components) == 2


def test_init():
    # Simple properties
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = GenericAbsSightline(radec)
    # Test
    assert gensys.sl_type == 'Generic'




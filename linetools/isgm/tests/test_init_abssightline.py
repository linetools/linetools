# Module to run tests on generating AbsSightline

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import GenericAbsSightline
from .utils import lyman_comp, si2_comp


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


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




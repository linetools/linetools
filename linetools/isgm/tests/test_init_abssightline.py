# Module to run tests on generating AbsSightline

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import GenericAbsSightline
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def lyman_comp(radec):
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lya.attrib['N'] = 1e17 /  u.cm**2
    lya.attrib['coord'] = radec
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    lyb.attrib['coord'] = radec
    abscomp = AbsComponent.from_abslines([lya,lyb])
    return abscomp


def si2_comp(radec):
    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.attrib['coord'] = radec
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines)
    #
    return SiII_comp


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




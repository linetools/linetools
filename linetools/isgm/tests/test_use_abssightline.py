# Module to run tests on generating AbsSightline

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

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
    SiII_comp.logN = 15.
    SiII_comp.flag_N = 1
    #
    return SiII_comp


def make_gensl():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    gensl = GenericAbsSightline.from_components([abscomp, SiII_comp])
    return gensl


def test_build_table():
    # Init
    gensl = make_gensl()
    # Table
    tbl = gensl.build_table()
    # Test
    idx = tbl['Ion'] == 'SiII'
    assert tbl['flag_N'][idx] == 1


def test_to_dict():
    # Init
    gensl = make_gensl()
    # Dict
    gensl_dict = gensl.to_dict()
    assert gensl_dict['class'] == 'GenericAbsSightline'

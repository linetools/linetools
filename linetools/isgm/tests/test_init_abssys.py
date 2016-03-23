# Module to run tests on generating AbsSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import GenericAbsSystem, LymanAbsSystem, AbsSystem
from linetools.spectralline import AbsLine

import pdb

class DLASystem(AbsSystem):
    def __init__(self, radec, zabs, vlim, NHI, **kwargs):
        AbsSystem.__init__(self, 'DLA', radec, zabs, vlim, NHI=NHI, **kwargs)
        pass

def test_init():
    # Simple properties
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = GenericAbsSystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=16.)
    # Test
    assert gensys.abs_type == 'Generic'
    np.testing.assert_allclose(gensys.zabs,1.244)


def test_init_strradec():
    # Simple properties
    gensys = GenericAbsSystem(('01:32:21', '+22:15:53.3'), 1.244, [-500,500]*u.km/u.s, NHI=16.)
    # Test
    assert gensys.abs_type == 'Generic'
    np.testing.assert_allclose(gensys.coord.ra.value, 23.087499999999995)
    #
    gensys = GenericAbsSystem('013221+221553.3', 1.244, [-500,500]*u.km/u.s, NHI=16.)
    np.testing.assert_allclose(gensys.coord.ra.value, 23.087499999999995)


def test_one_component():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # Instantiate
    HIsys = LymanAbsSystem.from_components([abscomp])
    # Test
    assert HIsys.abs_type == 'HILyman'
    assert len(HIsys._components) == 1
    assert HIsys._components[0].Zion[0] == 1
    assert HIsys._components[0].Zion[1] == 1


def test_multi_components():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines)
    SiII_comp.coord = radec
    # Instantiate
    LLSsys = GenericAbsSystem.from_components([abscomp,SiII_comp])
    # Test
    assert len(LLSsys._components) == 2


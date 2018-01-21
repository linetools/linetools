# Module to run tests on generating AbsSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import warnings

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import GenericAbsSystem, LymanAbsSystem
from linetools.isgm import utils as ltiu
from linetools.spectralline import AbsLine
from linetools.isgm.tests.utils import lyman_comp, si2_comp, oi_comp

import pdb

radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_from_json():
    # Tests from_dict too
    HIsys = LymanAbsSystem.from_json(data_path('HILya_abssys.json'))
    np.testing.assert_allclose(HIsys.zabs, 2.92939)
    assert len(HIsys._components[0].sig_logN) == 2
    # Tests ordering
    gensys = GenericAbsSystem.from_json(data_path('generic_abssys.json'))
    np.testing.assert_allclose(gensys._components[0].zcomp, 2.92939)


def test_write_json():
    HIsys = LymanAbsSystem.from_json(data_path('HILya_abssys.json'))
    HIsys.write_json()


def test_init():
    # Simple properties
    gensys = GenericAbsSystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=16.)
    # Test
    assert gensys.abs_type == 'Generic'
    assert gensys.flag_NHI == 1
    np.testing.assert_allclose(gensys.zabs,1.244)
    np.testing.assert_allclose(gensys.vlim.value, [-500,500])
    #


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
    # Instantiate
    abscomp = lyman_comp(radec)
    HIsys = LymanAbsSystem.from_components([abscomp])
    # Test
    assert HIsys.abs_type == 'HILyman'
    assert len(HIsys._components) == 1
    assert HIsys._components[0].Zion[0] == 1
    assert HIsys._components[0].Zion[1] == 1


def test_multi_components():
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    # Instantiate
    LLSsys = GenericAbsSystem.from_components([abscomp,SiII_comp])
    # Test
    np.testing.assert_allclose(LLSsys.NHI, 17.0)
    assert len(LLSsys._components) == 2


def test_add_component():
    abscomp = lyman_comp(radec)
    # Instantiate
    abssys = GenericAbsSystem.from_components([abscomp])
    # New component
    oicomp = oi_comp(radec, vlim=[-300.,300.]*u.km/u.s, z=abscomp.zcomp)
    # Standard
    assert abssys.add_component(oicomp)
    # Fail
    abssys = GenericAbsSystem.from_components([abscomp])
    oicomp2 = oi_comp(radec, vlim=[-400.,300.]*u.km/u.s, z=abscomp.zcomp)
    warnings.filterwarnings("ignore")
    assert not abssys.add_component(oicomp2)
    # Overlap
    assert abssys.add_component(oicomp2, overlap_only=True)


def test_build_systems_from_comp():
    # Build AbsSystems from components
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    #
    radec2 = SkyCoord(ra=223.1143*u.deg, dec=-12.4321*u.deg)
    abscomp2 = lyman_comp(radec2)
    #
    abs_systems = ltiu.build_systems_from_components([abscomp,SiII_comp,abscomp2])
    assert len(abs_systems) == 2


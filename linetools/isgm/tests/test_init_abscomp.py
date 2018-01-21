# Module to run tests on generating AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
import numpy as np
import pytest

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
import os

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_init():
    # Simple properties
    abscomp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)
    # Test
    assert abscomp.Zion[0] == 14
    np.testing.assert_allclose(abscomp.zcomp,1.0)
    print(abscomp)
    # Fine structure
    comp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s, Ej=0.1/u.cm) # need stars!
    assert comp.name.count('*') == 1


def test_init_failures():
    with pytest.raises(IOError):
        AbsComponent.from_abslines('blah')
    with pytest.raises(IOError):
        AbsComponent.from_abslines(['blah'])
    with pytest.raises(IOError):
        AbsComponent.from_component('blah')

    # Inconsistent abslines with median
    lya = AbsLine(1215.670*u.AA,z=2.92939)
    lya.attrib['N'] = 1e12 / u.cm**2
    lya.attrib['sig_N'] = [1e11]*2 / u.cm**2
    lya.attrib['b'] = 30 * u.km/u.s
    lyb = AbsLine('HI 1025',z=2.92939)
    lyb.attrib['N'] = 3e12 / u.cm**2
    lyb.attrib['sig_N'] = [3e11]*2 / u.cm**2
    lyb.attrib['b'] = 30 * u.km/u.s
    with pytest.raises(ValueError):
        AbsComponent.from_abslines([lya,lyb], adopt_median=True, chk_meas=True)

def test_read_from_json():
    abscomp = AbsComponent.from_json(data_path('old_lya_component.json'))
    assert np.isclose(abscomp.logN, 17.)
    assert abscomp.Zion == (1,1)
    assert len(abscomp.sig_logN) == 2
    assert abscomp.sig_N.size == 2


def test_init_single_absline():
    # Single AbsLine
    lya = AbsLine(1215.670*u.AA,z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lya.attrib['N'] = 1e12 / u.cm**2
    lya.attrib['sig_N'] = [1e11]*2 / u.cm**2
    lya.attrib['flag_N'] = 1
    abscomp = AbsComponent.from_abslines([lya])
    # Test
    assert abscomp.Zion[0] == 1
    assert len(abscomp.sig_N) == 2
    assert np.isclose(abscomp.sig_logN[0], 0.04342945)
    assert isinstance(abscomp.sig_logN, np.ndarray)
    np.testing.assert_allclose(abscomp.zcomp,2.92939)


def test_copy():
    # Single AbsLine
    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)
    abscomp = AbsComponent.from_abslines([lya])
    # Copy
    abscomp2 = abscomp.copy()
    # Checks
    attrs = vars(abscomp).keys()
    attrs2 = vars(abscomp2).keys()
    for attr in attrs:
        assert attr in attrs2
    np.testing.assert_allclose(abscomp._abslines[0].z,
                               abscomp2._abslines[0].z)


def test_init_multi_absline():
    # AbsLine(s)
    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lyb = AbsLine(1025.7222*u.AA)
    lyb.setz(lya.z)
    lyb.limits.set([-300.,300.]*u.km/u.s)
    # Instantiate
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # Test
    assert len(abscomp._abslines) == 2
    np.testing.assert_allclose(abscomp.zcomp,2.92939)

    # With column densities
    lya.attrib['N'] = 1e12 / u.cm**2
    lya.attrib['sig_N'] = [1e11]*2 / u.cm**2
    lya.attrib['flag_N'] = 1
    lya.attrib['b'] = 30 * u.km/u.s
    lyb.attrib['N'] = 3e12 / u.cm**2
    lyb.attrib['sig_N'] = [2e11]*2 / u.cm**2
    lyb.attrib['flag_N'] = 1
    lyb.attrib['b'] = 30 * u.km/u.s
    abscomp = AbsComponent.from_abslines([lya, lyb])
    # Test
    assert abscomp.flag_N == 1
    assert abscomp.attrib['sig_logN'].size == 2
    assert np.isclose(abscomp.logN, 12.146128035678238)

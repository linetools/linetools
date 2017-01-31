# Module to run tests on generating AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
import numpy as np
import pytest

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine


def test_init():
    # Simple properties
    abscomp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)
    # Test
    assert abscomp.Zion[0] == 14
    np.testing.assert_allclose(abscomp.zcomp,1.0)
    print(abscomp)


def test_init_failures():
    with pytest.raises(IOError):
        AbsComponent.from_abslines('blah')
    with pytest.raises(IOError):
        AbsComponent.from_abslines(['blah'])
    with pytest.raises(IOError):
        AbsComponent.from_component('blah')
    with pytest.raises(IOError):
        AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s, Ej=0.1/u.cm) # need stars!


def test_init_single_absline():
    # Single AbsLine
    lya = AbsLine(1215.670*u.AA,z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)
    abscomp = AbsComponent.from_abslines([lya])
    # Test
    assert abscomp.Zion[0] == 1
    np.testing.assert_allclose(abscomp.zcomp,2.92939)
    print(abscomp)


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
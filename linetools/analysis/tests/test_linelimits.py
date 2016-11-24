# Tests of LineLimits class
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pytest

from astropy import units as u

from linetools.analysis.linelimits import LineLimits
from linetools.spectralline import AbsLine

def test_init():
    # Init
    zlim = (0.999, 1.001)
    llim = LineLimits(1215.67*u.AA, 1., zlim)
    # Test
    with pytest.raises(AttributeError):
        llim.zlim=3
    # AbsLine
    lya = AbsLine('HI 1215')
    z=1.
    llim = LineLimits.from_specline(lya, z, zlim)
    # Bad zlim
    with pytest.raises(IOError):
        llim = LineLimits(1215.67*u.AA, 1., (1.1,1.2), chk_z=True)
    # Null zlim
    llim = LineLimits(1215.67*u.AA, 1., (1.,1))
    assert llim.is_set() is False

def test_set():
    # Init
    zlim = (0.999, 1.001)
    llim = LineLimits(1215.67*u.AA, 1., zlim)
    # Set
    llim.set(zlim)
    np.testing.assert_allclose(llim.wvlim.value, [2430.12433, 2432.55567])
    llim.set([2430.,2433.]*u.AA)
    np.testing.assert_allclose(llim.vlim.value, [-165.27207033, 204.61374917])
    llim.set([-160., 200]*u.km/u.s)

def test_use():
    # Init
    llim = LineLimits(1215.67*u.AA, 1., (0.999, 1.001))
    # Use
    np.testing.assert_allclose(llim.zlim, (0.999, 1.001))
    np.testing.assert_allclose(llim.wvlim.value, [2430.12433, 2432.55567])
    np.testing.assert_allclose(llim.vlim.value, [-149.93370, 149.858755])
    assert llim.vlim.unit == u.km/u.s
    # Print
    print(llim)

def test_to_dict():
    # Init
    llim = LineLimits(1215.67*u.AA, 1., (0.999, 1.001))
    ldict = llim.to_dict()
    for key in ['vlim', 'wrest', 'wvlim', 'z', 'zlim']:
        assert key in ldict.keys()

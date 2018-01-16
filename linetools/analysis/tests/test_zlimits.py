# Tests of zLimits class
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pytest

from astropy import units as u

from linetools.analysis.zlimits import zLimits
from linetools.spectralline import AbsLine

def test_init():
    # Init
    zlim = (0.999, 1.001)
    llim = zLimits(1., zlim)
    # Init
    zlim = (0.999, 1.001)
    llim = zLimits(1., zlim, wrest=1215.67*u.AA)
    # Test
    with pytest.raises(AttributeError):
        llim.zlim=3
    # AbsLine
    lya = AbsLine('HI 1215')
    z=1.
    zlim = (0.999, 1.001)
    llim = zLimits.from_specline(lya, z, zlim)
    # Bad zlim
    with pytest.raises(IOError):
        llim = zLimits(1., (1.1,1.2), wrest=1215.67*u.AA, chk_z=True)
    # Null zlim
    llim = zLimits(1., (1.,1), wrest=1215.67*u.AA)
    assert llim.is_set() is False

def test_set():
    # Init
    zlim = (0.999, 1.001)
    llim = zLimits(1., zlim, wrest=1215.67*u.AA)
    # Set
    llim.set(zlim)
    np.testing.assert_allclose(llim.wvlim.value, [2430.12433, 2432.55567])
    llim.set([2430.,2433.]*u.AA)
    np.testing.assert_allclose(llim.vlim.value, [-165.27207033, 204.61374917])
    llim.set([-160., 200]*u.km/u.s)

def test_use():
    # Init
    llim = zLimits(1., (0.999, 1.001), wrest=1215.67*u.AA)
    # Use
    np.testing.assert_allclose(llim.zlim, (0.999, 1.001))
    np.testing.assert_allclose(llim.wvlim.value, [2430.12433, 2432.55567])
    np.testing.assert_allclose(llim.vlim.value, [-149.93370, 149.858755])
    np.testing.assert_allclose(llim.vmin.value, -149.93370)
    np.testing.assert_allclose(llim.vmax.value, 149.858755)
    assert llim.vlim.unit == u.km/u.s
    # Print
    print(llim)

def test_to_dict():
    # Init
    llim = zLimits(1., (0.999, 1.001), wrest=1215.67*u.AA)
    ldict = llim.to_dict()
    for key in ['vlim', 'wrest', 'wvlim', 'z', 'zlim']:
        assert key in ldict.keys()

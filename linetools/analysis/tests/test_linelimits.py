# Tests of LineLimits class
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest

from astropy import units as u

from linetools.analysis.linelimits import LineLimits

def test_init():
    # Init
    llim = LineLimits(1215.67*u.AA, 1., (0.999, 1.001))
    # Test
    with pytest.raises(AttributeError):
        llim.zlim=3

def test_use():
    # Init
    llim = LineLimits(1215.67*u.AA, 1., (0.999, 1.001))
    # Use
    np.testing.assert_allclose(llim.zlim, (0.999, 1.001))
    np.testing.assert_allclose(llim.wvlim.value, [2430.12433, 2432.55567])
    np.testing.assert_allclose(llim.vlim.value, [-149.896229, 149.896229])
    assert llim.vlim.unit == u.km/u.s

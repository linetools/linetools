# Tests of LineLimits class
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest

from astropy import units as u

from linetools.analysis.linelimits import LineLimits

def test_init():
    # Make fake spectrum
    llim = LineLimits(1215.67*u.AA, 1., (0.999, 1.001))
    # Test
    #np.testing.assert_allclose((N.value, sig_N.value),
    #                           (96652191688169.72, 194151305045168.12))
    #assert N.unit == u.cm**-2


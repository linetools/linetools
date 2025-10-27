from __future__ import print_function, absolute_import, division, unicode_literals

from ...spectra.io import readspec
from ..continuumfnd import contknots
import os
import numpy as np
import pytest

def test_find_continuum_knots():
    # Get the linetools package directory from this module's location
    d = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    spec = readspec(d + '/spectra/tests/files/spec_example_2.fits')  # XSpectrum1D.from_file(spfile)
    testknots, testknotpixs = contknots(spec, ewsnlim=5, showcont=False)
    knotscont = np.asarray(testknots)[:, 1]
    i = np.where(np.isnan(knotscont) == False)
    dy = abs(knotscont[i] - 1.)
    assert len(dy) > 5
    assert np.allclose(dy, np.repeat(0., len(dy)), rtol=0., atol=0.1)
    assert np.std(dy) < 0.02


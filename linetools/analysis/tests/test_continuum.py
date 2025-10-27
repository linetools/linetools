from __future__ import print_function, absolute_import, division, unicode_literals

from astropy.table import Table
from ...spectra.io import readspec
from ..continuum import find_continuum
import os
import numpy as np
import pytest

def test_find_continuum():
    # Get the linetools package directory from this module's location
    d = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    spec = readspec(d + '/spectra/tests/files/q0002m422.txt.gz', masking='none')
    co, pts = find_continuum(spec, redshift=2.76, divmult=3.5,
                       forest_divmult=3, kind='QSO')
    assert np.allclose(co[:3],
                       [19624.02210835, 19615.23362468, 19606.44507182])
    assert np.allclose(co[-3:],[384.69366713, 384.69095479,   384.68824243])
    assert np.allclose(co[50000:50003],
                       [1575.88179631, 1575.5837134, 1575.28703315])

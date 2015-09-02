from __future__ import print_function, absolute_import, division, unicode_literals

from astropy.table import Table
from ...spectra.io import readspec
from ..continuum import find_continuum
import imp
import numpy as np

def test_find_continuum():
    d = imp.find_module('linetools')[1]
    spec = readspec(d + '/spectra/tests/files/q0002m422.txt.gz')
    co, pts = find_continuum(spec, redshift=2.76, divmult=3.5,
                       forest_divmult=3, kind='QSO')
    assert np.allclose(co[:3], [6810.73400022,  6807.63940372,  6804.54478286])
    assert np.allclose(co[-3:],[384.69366713, 384.69095479,   384.68824243])
    assert np.allclose(co[50000:50003],
                       [1575.88179631, 1575.5837134, 1575.28703315])

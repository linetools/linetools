from __future__ import print_function, absolute_import, division, unicode_literals

from ..interp import interp_Akima
import numpy as np

def test_interp_Akima():
    x = np.sort(np.random.random(10) * 10)
    y = np.random.normal(0.0, 0.1, size=len(x))
    assert np.allclose(y, interp_Akima(x, x, y))

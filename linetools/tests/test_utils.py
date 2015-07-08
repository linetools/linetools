from __future__ import print_function, absolute_import, division, unicode_literals

from ..utils import between
import numpy as np

def test_between():
    x = [1,2,3,4,5]
    c0 = between(x, 2.5, 3.5)
    assert np.all(c0 == np.array([False, False, True, False, False]))

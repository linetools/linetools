from __future__ import print_function, absolute_import, division, unicode_literals

from ..utils import between, radec_to_coord
import numpy as np
import pdb

def test_between():
    x = [1,2,3,4,5]
    c0 = between(x, 2.5, 3.5)
    assert np.all(c0 == np.array([False, False, True, False, False]))

def test_radeccoord():
    for radec in ['J124511+144523', '124511+144523',
                  'J12:45:11+14:45:23', ('12:45:11','+14:45:23')]:
        coord = radec_to_coord(radec)
        # Test
        np.testing.assert_allclose(coord.ra.value, 191.2958333333333)

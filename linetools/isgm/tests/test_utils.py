# Module to run tests utils on isigm

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

from linetools.isgm import utils as ltiu
import astropy.units as u
import numpy as np

import imp
lt_path = imp.find_module('linetools')[1]


def test_overlapping_chunks():
    chunk1 = (1,2,3,4)
    chunk2 = [3,4,5,6]
    t = ltiu.overlapping_chunks(chunk1, chunk2)
    assert t
    t = ltiu.overlapping_chunks(chunk2, chunk1)
    assert t

    chunk2 = np.array([5,7])
    f = ltiu.overlapping_chunks(chunk2, chunk1)
    assert ~f
    f = ltiu.overlapping_chunks(chunk2*u.AA, chunk1*u.AA)
    assert ~f

    # Wrong format
    try:
        f = ltiu.overlapping_chunks(chunk2*u.AA, chunk1)
    except ValueError:
        pass
    try:
        f = ltiu.overlapping_chunks(chunk2*u.AA, chunk1*u.K)
    except ValueError:
        pass

    # not sorted
    try:
        f = ltiu.overlapping_chunks((1,2,3,4), (4,5,3,6))
    except ValueError:
        pass




# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy.coordinates import SkyCoord

from linetools.isgm.io import read_joebvp

#import pdb
#pdb.set_trace()
# Set of Input lines

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_read_joebvp():
    vp_file = data_path('group_9.VP')
    # Load
    icoord = SkyCoord(ra=12., dec=-12, unit='deg')
    comps = read_joebvp(vp_file, icoord)
    # Test
    assert isinstance(comps, list)
    assert len(comps) == 2
    assert len(comps[0]._abslines) == 2



# Module to run tests on SolarAbund class

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pdb
import pytest
#from astropy import units as u
import numpy as np

from linetools.abund import solar
#import pdb
#pdb.set_trace()
# Set of Input lines

# Simple init
def test_init():
    sol = solar.SolarAbund()
    assert sol.ref == 'Asplund2009'
    print(sol)
    # Failures
    pytest.raises(TypeError, solar.SolarAbund, 1)
    pytest.raises(ValueError, solar.SolarAbund, ref='dum')


def test_elm():
    sol = solar.SolarAbund()
    np.testing.assert_allclose(sol['C'],8.43)
    np.testing.assert_allclose(sol[['C', 'Fe']],np.array([8.43,7.45]))


def test_Z():
    sol = solar.SolarAbund()
    np.testing.assert_allclose(sol[6],8.43)


def test_ratio():
    sol = solar.SolarAbund()
    np.testing.assert_allclose(sol.get_ratio('C/Fe'),0.98)


def test_getitem_errors():
    sol = solar.SolarAbund()
    with pytest.raises(ValueError):
        a = sol[-1]
    with pytest.raises(ValueError):
        a = sol['bad_name']
    with pytest.raises(IndexError):
        a = sol[None]  # bad input




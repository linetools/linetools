# Module to run tests on SolarAbund class

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import os
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
	# 
	assert sol.ref == 'Asplund2009'

def test_elm():
	sol = solar.SolarAbund()
	np.testing.assert_allclose(sol['C'],8.43)

def test_Z():
	sol = solar.SolarAbund()
	np.testing.assert_allclose(sol[6],8.43)

def test_ratio():
	sol = solar.SolarAbund()
	np.testing.assert_allclose(sol['C/Fe'],0.98)



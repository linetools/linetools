# Module to run tests on Generating a LineList
#   Also tests some simple functionality

# TEST_UNICODE_LITERALS

import os, pdb
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np

from linetools.lists import parse

# Morton 2003 ASCII file
def test_morton03():
	m03 = parse.parse_morton03(orig=True)
	# 
	np.testing.assert_allclose(m03['wrest'][5], 1025.7218, rtol=1e-7)

	assert m03['wrest'].unit == u.Angstrom

# Morton 2000 ASCII file
def test_morton00():
	m00 = parse.parse_morton00(orig=True)

	np.testing.assert_allclose(m00['wrest'][5], 2593.3093, rtol=1e-7)

	assert m00['wrest'].unit == u.AA

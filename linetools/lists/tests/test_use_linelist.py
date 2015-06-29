# Module to run tests on Generating a LineList
#   Also tests some simple functionality

# TEST_UNICODE_LITERALS

import os
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np

from linetools.lists.linelist import LineList

#import pdb
#pdb.set_trace()

# ISM LineList
def test_closest():
	ism = LineList('ISM')
	ism.closest=True
	# 
	line = ism[1250.584*u.AA]
	np.testing.assert_allclose(line['wrest'], 1250.578*u.AA, rtol=1e-7)


# Module to run tests on Generating a LineList
#   Also tests some simple functionality

# TEST_UNICODE_LITERALS

import os, pdb
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np

from linetools.lists.linelist import LineList

#import pdb
#pdb.set_trace()

# ISM LineList
def test_lines_from_ion():
	ism = LineList('ISM')
	# 
	lines = ism[(6,2)]
	assert (1334.5323*u.AA in lines['wrest'])

def test_closest():
	ism = LineList('ISM')
	ism.closest=True
	# 
	line = ism[1250.584*u.AA]
	np.testing.assert_allclose(line['wrest'], 1250.578*u.AA, rtol=1e-7)

def test_all_transition_names_from_line():
	ism = LineList('ISM')
	line = 'OVI'
	names = ism.all_transition_names_from_line(line)
	assert len(names == 2), 'Something is wrong in all_transition_names_from_line()'
	print(names)

# Module to run tests on ion code

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from linetools.abund import ions

#import pdb
#pdb.set_trace()
# Set of Input lines

def test_ion_to_name():
	# Normal
	ionnm = ions.ion_name((6,2)) 
	assert ionnm == 'CII'
	# Latex
	ionnm = ions.ion_name((6,2),flg=1) 
	assert ionnm == '{\\rm C}^{+}'

def test_name_to_ion():
	Zion = ions.name_ion('Si II')
	assert Zion == (14,2)

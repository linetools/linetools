# Module to run tests on LineList utilities

from __future__ import (print_function, absolute_import, division,
unicode_literals)

# TEST_UNICODE_LITERALS

import pytest

from astropy import units as u

from linetools.lists.linelist import LineList
from .. import utils

#import pdb
#pdb.set_trace()

# ISM LineList
def test_dict_to_table():
    #
    dct = dict(wrest=1215.67*u.AA, f=0.12, A=1e-6/u.s)
    tab = utils.from_dict_to_table(dct)
    # Check
    assert tab['wrest'].unit == u.AA

def test_table_to_dict():
    ism = LineList('ISM')
    #
    dct = utils.from_table_to_dict(ism._data[0:1])
    # Check units
    assert dct['wrest'].unit == u.AA


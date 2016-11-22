# Module to run tests on generating AbsSightline

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
from .utils import make_gensl


def test_build_table():
    # Init
    gensl = make_gensl()
    # Table
    tbl = gensl.build_table()
    # Test
    idx = tbl['Z'] == 14
    assert tbl['flag_N'][idx] == 1


def test_to_dict():
    # Init
    gensl = make_gensl()
    # Dict
    gensl_dict = gensl.to_dict()
    assert gensl_dict['class'] == 'GenericAbsSightline'

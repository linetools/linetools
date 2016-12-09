# Tests of linetools.utils

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
import os

import linetools.io as lio


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_emline_from_alis():
    lio.emlines_from_alis_output(data_path('spec1d_J0018p2345_KASTb_coadd.mod.out'))


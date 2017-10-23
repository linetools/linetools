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
    J0018_emlines = lio.emlines_from_alis_output(data_path('spec1d_J0018p2345_KASTb_coadd.mod.out'))

    assert len(J0018_emlines) == 8
    np.testing.assert_allclose(J0018_emlines[0].z, 0.01540425)
    np.testing.assert_allclose(J0018_emlines[1].wrest.value, 4341.684) #Hgamma


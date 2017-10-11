# Module to run tests on generating AbsSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.emsystem import EmSystem

coord = "00:18:59.32 +23:45:40.32"
radec = SkyCoord(coord, unit=(u.hourangle, u.deg))

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():
    emsys = EmSystem(radec, 0.0154)
    assert np.isclose(emsys.zem, 0.0154)


def test_from_alis():
    # Tests from_dict too
    alis_file = data_path('spec1d_J0018p2345_KASTb_coadd.mod.out')
    emsys = EmSystem.from_alis(alis_file, 'J001859.32+234540.32')
    #
    assert np.isclose(emsys.zem, 0.01540425)


def test_add_lines_from_alis():
    # Tests from_dict too
    alisb_file = data_path('spec1d_J0018p2345_KASTb_coadd.mod.out')
    emsys = EmSystem.from_alis(alisb_file, 'J001859.32+234540.32')
    #
    alisr_file = data_path('spec1d_J0018p2345_KASTr_coadd.mod.out')
    emsys.add_emlines_from_alis(alisr_file, chk_z=False)
    assert len(emsys._emlines) == 13

def test_emsys_to_dict():
    emsys = EmSystem(radec, 0.0154)
    todict = emsys.to_dict()

    assert len(todict) == 14
    np.isclose(todict['RA'],4.747166666666666)
    np.isclose(todict['DEC'], 23.7612)
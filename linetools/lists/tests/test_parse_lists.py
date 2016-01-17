# Module to run tests on Generating a LineList
#   Also tests some simple functionality

# TEST_UNICODE_LITERALS

import pdb
import pytest
from astropy import units as u
import numpy as np

from linetools.lists import parse

# Morton 2003 ASCII file

@pytest.mark.skipif("sys.version_info >= (3,0)")
def test_morton03():
    m03 = parse.parse_morton03(orig=True)
    #
    np.testing.assert_allclose(m03['wrest'][5], 930.7482, rtol=1e-7)

    assert m03['wrest'].unit == u.Angstrom
    #
    m03 = parse.mktab_morton03()
    #m03 = parse.mktab_morton03(do_this=True, outfil='tmp.fits')
    m03 = parse.mktab_morton03(do_this=True, fits=False, outfil='tmp.vo')

# Morton 2000 ASCII file
def test_morton00():
    m00 = parse.parse_morton00(orig=True)

    np.testing.assert_allclose(m00['wrest'][5], 2593.3093, rtol=1e-7)

    assert m00['wrest'].unit == u.AA
    #
    m00 = parse.mktab_morton00()
    #m00 = parse.mktab_morton00(do_this=True, outfil='tmp.fits')

def test_verner96():
    v96 = parse.parse_verner96(orig=True)

def test_galaxy_lines():
    glx = parse.grab_galaxy_linelists()

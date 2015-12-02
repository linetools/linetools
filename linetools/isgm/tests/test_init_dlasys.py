# Module to run tests on initializing DLA System

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from linetools.isgm.dla import DLASystem, DLASurvey

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_simple_dla_init():
	# Init 
    dla = DLASystem((0.*u.deg, 0.*u.deg), 2.5, None, NHI=20.55)
    #
    np.testing.assert_allclose(dla.vlim[0].value,-500.)
    np.testing.assert_allclose(dla.NHI, 20.55)

def test_dat_init():
    # JXP .dat files
    if os.getenv('DLA') is None:
        assert True
        return
    # Read
    datfil = 'Data/PH957.z2309.dat'
    dla = DLASystem.from_datfile(datfil, tree=os.environ.get('DLA'))
    #
    np.testing.assert_allclose(dla.NHI, 21.37)
    np.testing.assert_allclose(dla.zabs, 2.309)

def test_parse_ion():
    # JXP .ion file
    if os.getenv('DLA') is None:
        assert True
        return
    # Read
    datfil = 'Data/PH957.z2309.dat'
    dla = DLASystem.from_datfile(datfil, tree=os.environ.get('DLA'))
    #
    dla.get_ions(use_clmfile=True)
    assert len(dla._ionN) == 13

def test_default_dla_sample():
    if os.getenv('DLA') is None:
        assert True
        return
    # Load
    dlas = DLASurvey.default_sample()
    assert len(dlas._abs_sys) == 100

def test_default_dla_sample_with_ions():
    if os.getenv('DLA') is None:
        assert True
        return
    # Load
    dlas = DLASurvey.default_sample()
    dlas.fill_ions(use_clmfile=True)
    CIV_clms = dlas.ions((6,4))
    gdCIV = np.where(CIV_clms['flag_N']>0)[0]
    assert len(gdCIV) == 74


# Module to run tests on initializing DLA System

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from linetools.isgm.lls import DLASystem

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
    if os.getenv('LLSTREE') is None:
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
    if os.getenv('LLSTREE') is None:
        assert True
        return
    # Read
    datfil = 'Data/UM184.z2929.dat'
    lls = LLSSystem.from_datfile(datfil, tree=os.getenv('LLSTREE'))
    #
    dla.get_ions(use_clmfile=True)
    assert len(lls._ionN) == 13
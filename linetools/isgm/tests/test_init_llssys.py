# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from xastropy.igm.abs_sys.lls_utils import LLSSystem

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_simple_init():
	# Init 
    lls = LLSSystem((0.*u.deg, 0.*u.deg), 2.0, None, NHI=17.9)
    #
    np.testing.assert_allclose(lls.NHI, 17.9)
    np.testing.assert_allclose(lls.tau_LL, 5.03537353413629)

def test_dat_init():
    # JXP .dat files
    if os.getenv('LLSTREE') is None:
        assert True
        return
    # Read
    datfil = 'Data/UM184.z2929.dat'
    lls = LLSSystem.from_datfile(datfil, tree=os.getenv('LLSTREE'))
    #    
    np.testing.assert_allclose(lls.zabs, 2.93012)

def test_parse_ion():
    # JXP .ion file
    if os.getenv('LLSTREE') is None:
        assert True
        return
    # Read
    datfil = 'Data/UM184.z2929.dat'
    lls = LLSSystem.from_datfile(datfil, tree=os.getenv('LLSTREE'))
    #    
    lls.get_ions(use_clmfile=True)
    assert len(lls._ionclms) == 13
# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import glob, os, imp, pdb
import pytest
#from astropy import units as u
#from astropy.coordinates import SkyCoord

from xastropy.igm.abs_sys.lls_utils import LLSSurvey

xa_path = imp.find_module('xastropy')[1]

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_read_hdlls_dr1():
    # Read Summary
    summ_fil = glob.glob(xa_path+"/data/LLS/HD-LLS_DR1.fits")
    if len(summ_fil) > 0:
        summ_fil = summ_fil[0]
    else:
        assert True
        return
    lls = LLSSurvey.from_sfits(summ_fil)
    assert lls.nsys == 157

    # Read ions
    ions_fil = glob.glob(xa_path+"/data/LLS/HD-LLS_ions.json")
    if len(ions_fil) > 0:
        ions_fil = ions_fil[0]
    else:
        assert True
        return
    lls.fill_ions(jfile=ions_fil)
    CII_clms = lls.ions((6,2))
    gdCII = np.where(CII_clms['flag_N']>0)[0]
    assert len(gdCII) == 103

def test_dat_list():
    '''JXP format :: Likely to be Deprecated
    '''
    # LLS Survey
    if os.getenv('LLSTREE') is None:
        assert True
        return
    # Load
    lls = LLSSurvey.from_flist('Lists/lls_metals.lst', tree=os.getenv('LLSTREE'))
    # tests
    np.testing.assert_allclose(lls.NHI[0], 19.25)
    assert lls.nsys == 165


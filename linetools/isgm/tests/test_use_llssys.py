# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
from astropy import units as u

from linetools.spectra.xspectrum1d import XSpectrum1D

from xastropy.igm.abs_sys.lls_utils import LLSSystem

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_fluxmodel():
	# Init
    lls = LLSSystem((0.*u.deg, 0.*u.deg), 2.5, None, NHI=17.9)
    # Fill LLS lines
    lls.fill_lls_lines()
    # Generate a spectrum
    wave = np.arange(3000., 6500)
    npix = len(wave)
    spec = XSpectrum1D.from_tuple((wave*u.AA,np.ones(npix)))
    # Model
    model = lls.flux_model(spec)
    np.testing.assert_allclose(model.flux[100].value,0.009424664763760516)


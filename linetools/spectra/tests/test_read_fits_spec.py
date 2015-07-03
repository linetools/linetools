# Module to run tests on spectra.io
import os
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np
from astropy.io import fits

from linetools.spectra import io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

# Separate flux/error files
def test_sep_files():
    spec = io.readspec(data_path('UM184_nF.fits'))
    idl = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])
    np.testing.assert_allclose(spec.dispersion, idl['wave'])
    np.testing.assert_allclose(spec.sig, idl['sig'], atol=1e-3, rtol=0)

    assert spec.dispersion.unit == u.Unit('AA')

def test_setwave():
    hd = fits.getheader(data_path('UM184_nF.fits'))
    wave = io.setwave(hd)
    np.testing.assert_allclose(wave[0], 3042.34515916)
    hd['CRPIX1'] = 10
    wave = io.setwave(hd)
    np.testing.assert_allclose(wave[0], 3040.33648468)


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

# Rel vel
def test_addnoise():
    spec = io.readspec(data_path('UM184_nF.fits'))
    #
    spec.add_noise(seed=12)
    np.testing.assert_allclose(spec.flux[1000], 0.44806158542633057)

    # With S/N input
    spec.add_noise(seed=19,s2n=10.)
    np.testing.assert_allclose(spec.flux[1000], 0.24104823059199412)

# Box car smooth
def test_box_smooth():
    spec = io.readspec(data_path('UM184_nF.fits'))

    # Smooth
    newspec3 = spec.box_smooth(3)
    np.testing.assert_allclose(newspec3.flux[4000], 0.9967582821846008)
    assert newspec3.flux.unit == u.dimensionless_unscaled

    newspec5 = spec.box_smooth(5)
    np.testing.assert_allclose(newspec5.flux[3000], 1.086308240890503)

# Gaussian smooth
def test_gauss_smooth():
    spec = io.readspec(data_path('UM184_nF.fits'))

    # Smooth
    smth_spec = spec.gauss_smooth(4.)
    # Test
    np.testing.assert_allclose(smth_spec.flux[3000].value, 0.8288110494613)
    assert smth_spec.flux.unit == spec.flux.unit

# Rebin
def test_rebin():
    spec = io.readspec(data_path('UM184_nF.fits'))
    # Rebin
    new_wv = np.arange(3000., 9000., 5) * u.AA
    newspec = spec.rebin(new_wv)
    # Test

    np.testing.assert_allclose(newspec.flux[1000], 0.9999280967617779)
    assert newspec.flux.unit == u.dimensionless_unscaled

# Rel vel
def test_relvel():
    spec = io.readspec(data_path('UM184_nF.fits'))

    # Velocity
    velo = spec.relative_vel(5000.*u.AA)
    # Test
    np.testing.assert_allclose(velo[6600].value, -3716.441360213781)
    assert velo.unit == (u.km/u.s)

# Write FITS
def test_write_fits():
    spec = io.readspec(data_path('UM184_nF.fits'))

    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_fits(data_path('tmp.fits'))
    spec2 = io.readspec(data_path('tmp.fits'))
    # check a round trip works
    np.testing.assert_allclose(spec.dispersion, spec2.dispersion)

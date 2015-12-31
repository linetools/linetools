# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
import pdb
from astropy import units as u
import numpy as np
from astropy.io import fits, ascii
from astropy.table import QTable, Table, Column

from linetools.spectra import io
from linetools.spectra.xspectrum1d import XSpectrum1D

@pytest.fixture
def spec():
    return io.readspec(data_path('UM184_nF.fits'))


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_rebin(spec):
    # Rebin
    new_wv = np.arange(3000., 9000., 5) * u.AA
    newspec = spec.rebin(new_wv)
    # Test
    np.testing.assert_allclose(newspec.flux[1000], 0.9999280967617779)
    assert newspec.flux.unit == u.dimensionless_unscaled
    # With sigma
    newspec = spec.rebin(new_wv, do_sig=True)
    pdb.set_trace()
    np.testing.assert_allclose(newspec.flux[1000], 0.9999280967617779)

def test_addnoise(spec):
    #
    spec.add_noise(seed=12)
    np.testing.assert_allclose(spec.flux[1000], 0.44806158542633057)

    # With S/N input
    spec.add_noise(seed=19,s2n=10.)
    np.testing.assert_allclose(spec.flux[1000], 0.24104823059199412)


def test_box_smooth(spec):

    # Smooth
    newspec3 = spec.box_smooth(3)
    np.testing.assert_allclose(newspec3.flux[4000], 0.9967582821846008)
    assert newspec3.flux.unit == u.dimensionless_unscaled

    newspec5 = spec.box_smooth(5)
    np.testing.assert_allclose(newspec5.flux[3000], 1.086308240890503)


def test_gauss_smooth(spec):

    # Smooth
    smth_spec = spec.gauss_smooth(4.)
    # Test
    np.testing.assert_allclose(smth_spec.flux[3000].value, 0.82889723777)
    assert smth_spec.flux.unit == spec.flux.unit




def test_relvel(spec):

    # Velocity
    velo = spec.relative_vel(5000.*u.AA)
    # Test
    np.testing.assert_allclose(velo[6600].value, -3716.441360213781)
    assert velo.unit == (u.km/u.s)


def test_print_repr(spec):
    print(repr(spec))
    print(spec)


def test_write_ascii(spec):
    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_ascii(data_path('tmp.ascii'))
    # 
    spec2 = io.readspec(data_path('tmp.ascii'))
    # check a round trip works
    np.testing.assert_allclose(spec.dispersion, spec2.dispersion)


def test_write_fits(spec):
    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_fits(data_path('tmp.fits'))
    spec2 = io.readspec(data_path('tmp.fits'))
    # check a round trip works
    np.testing.assert_allclose(spec.dispersion, spec2.dispersion)


def test_readwrite_without_sig():
    sp = XSpectrum1D.from_tuple(([5,6,7], np.ones(3)))
    sp.write_to_fits(data_path('tmp.fits'))
    sp1 = io.readspec(data_path('tmp.fits'))
    np.testing.assert_allclose(sp1.dispersion.value, sp.dispersion.value)
    np.testing.assert_allclose(sp1.flux.value, sp.flux.value)


def test_readwrite_metadata(spec):
    d = {'a':1, 'b':'abc', 'c':3.2, 'd':np.array([1,2,3]),
         'e':dict(a=1,b=2)}
    spec.meta.update(d)
    spec.write_to_fits(data_path('tmp.fits'))
    spec2 = io.readspec(data_path('tmp.fits'))
    assert spec2.meta['a'] == d['a']
    assert spec2.meta['b'] == d['b']
    np.testing.assert_allclose(spec2.meta['c'], d['c'])
    np.testing.assert_allclose(spec2.meta['d'], d['d'])
    assert spec2.meta['e'] == d['e']


def test_copy(spec):
    spec2 = spec.copy()
    assert spec.wavelength[0] == spec2.wavelength[0]
    assert spec.flux[-1] == spec2.flux[-1]


def test_continuum_utils(spec):
    # define continuum in a non-interactive way...
    n = len(spec.wavelength.value)
    contpoints = [(spec.wavelength.value[i], 1.) for i in range(n)[::int(n/100)]]
    spec.meta['contpoints'] = contpoints
    xy = np.array(spec.meta['contpoints'])
    xy = xy.transpose()
    x, y = xy[0], xy[1]
    # test interpolate
    spec.co = spec._interp_continuum(x, y)
    np.testing.assert_allclose(spec.co, 1.)
    co_old = spec.co
    # test perturb
    spec.perturb_continuum(rel_var=0.05, seed=2)
    assert all(co_old != spec.co)
    np.testing.assert_allclose(np.max(spec.co), 1.2319426067564621)
    np.testing.assert_allclose(np.min(spec.co), 0.86589518482815)
    #test reset
    spec.reset_continuum()
    np.testing.assert_allclose(spec.co, 1.)

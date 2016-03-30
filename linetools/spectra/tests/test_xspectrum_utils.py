# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
import pdb
from astropy import units as u
import numpy as np

from linetools.spectra import io
from linetools.spectra.xspectrum1d import XSpectrum1D

@pytest.fixture
def spec():
    return io.readspec(data_path('UM184_nF.fits'))

@pytest.fixture
def spec2():
    return io.readspec(data_path('PH957_f.fits'))

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_rebin(spec):
    # Rebin
    new_wv = np.arange(3000., 9000., 5) * u.AA
    newspec = spec.rebin(new_wv)
    # Test
    np.testing.assert_allclose(newspec.flux[1000], 0.9999280967617779)
    assert newspec.flux.unit is u.dimensionless_unscaled
    # With sigma
    newspec = spec.rebin(new_wv, do_sig=True)
    """
    i1 = np.argmin(np.abs(spec.wavelength-5000.*u.AA))
    s2n_1 = spec.flux[i1] / spec.sig[i1]
    i2 = np.argmin(np.abs(newspec.wavelength-5000.*u.AA))
    s2n_2 = newspec.flux[i2] / newspec.sig[i2]
    """
    np.testing.assert_allclose(newspec.sig[1000].value, 0.01432187, rtol=1e-5)

def test_addnoise(spec):
    #
    newspec = spec.add_noise(seed=12)
    np.testing.assert_allclose(newspec.flux[1000].value, 0.4480615, rtol=1e-5)

    # With S/N input
    newspec2 = spec.add_noise(seed=19,s2n=10.)
    np.testing.assert_allclose(newspec2.flux[1000].value, -0.3028939, rtol=1e-5)


def test_box_smooth(spec):

    # Smooth
    newspec3 = spec.box_smooth(3)
    np.testing.assert_allclose(newspec3.flux[4000], 0.9967582821846008)
    assert newspec3.flux.unit == u.dimensionless_unscaled

    newspec5 = spec.box_smooth(5)
    np.testing.assert_allclose(newspec5.flux[3000], 1.086308240890503)
    # Preserve
    newspec5p = spec.box_smooth(5, preserve=True)


def test_gauss_smooth(spec):

    # Smooth
    smth_spec = spec.gauss_smooth(4.)
    # Test
    np.testing.assert_allclose(smth_spec.flux[3000].value, 0.82889723777)
    assert smth_spec.flux.unit == spec.flux.unit


def test_print_repr(spec):
    print(repr(spec))
    print(spec)


def test_rebin(spec):
    # Rebin
    new_wv = np.arange(3000., 9000., 5) * u.AA
    newspec = spec.rebin(new_wv, do_sig=True)
    # Test

    np.testing.assert_allclose(newspec.flux[1000], 0.9999280967617779)
    assert newspec.flux.unit == u.dimensionless_unscaled


def test_relvel(spec):

    # Velocity
    velo = spec.relative_vel(5000.*u.AA)
    # Test
    np.testing.assert_allclose(velo[6600].value, -3716.441360213781)
    assert velo.unit == (u.km/u.s)


def test_splice(spec, spec2):
    spec3 = spec.splice(spec2)


def test_write_ascii(spec):
    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_ascii(data_path('tmp.ascii'))
    #
    spec2 = io.readspec(data_path('tmp.ascii'))
    # check a round trip works
    np.testing.assert_allclose(spec.wavelength, spec2.wavelength)


def test_write_fits(spec, spec2):
    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_fits(data_path('tmp.fits'))
    specin = io.readspec(data_path('tmp.fits'))
    # check a round trip works
    np.testing.assert_allclose(spec.wavelength, specin.wavelength)
    # ESI
    spec2.write_to_fits(data_path('tmp2.fits'))
    specin2 = io.readspec(data_path('tmp2.fits'))
    # check a round trip works
    np.testing.assert_allclose(spec2.wavelength, specin2.wavelength)


def test_readwrite_without_sig():
    sp = XSpectrum1D.from_tuple((np.array([5,6,7]), np.ones(3)))
    sp.write_to_fits(data_path('tmp.fits'))
    sp1 = io.readspec(data_path('tmp.fits'))
    np.testing.assert_allclose(sp1.wavelength.value, sp.wavelength.value)
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

def test_plot(spec):
    spec.plot(show=False)

def test_continuum_utils(spec):
    # define continuum in a non-interactive way...
    n = len(spec.wavelength.value)
    contpoints = [(spec.wavelength.value[i], 1.) for i in range(n)[::int(n/100)]]
    spec.meta['contpoints'] = contpoints
    xy = np.array(spec.meta['contpoints'])
    xy = xy.transpose()
    x, y = xy[0], xy[1]
    # test interpolate
    spec.normalize(spec._interp_continuum(x, y, spec.wavelength.value))
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

    # test normalize/unnormalize
    flux_old = spec.flux
    spec.unnormalize()
    assert spec.normed is False
    np.testing.assert_allclose(spec.flux,flux_old)


def test_assignment(spec):
    temp = np.arange(1, len(spec.wavelength) + 1)
    spec.wavelength =  temp * u.m
    assert spec.wavelength[0] == temp[0] * u.m
    unit = u.erg / u.s / u.cm**2 / u.AA
    spec.flux = temp * unit
    assert spec.flux[0] == temp[0] * unit
    spec.sig = temp
    assert spec.sig[0] == temp[0] * unit
    spec.co = temp
    assert spec.co[0] == temp[0] * unit

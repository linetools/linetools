# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
from astropy import units as u
import numpy as np

from linetools.spectra import io
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import utils as ltsu

@pytest.fixture
def spec():
    return io.readspec(data_path('UM184_nF.fits'))


@pytest.fixture
def spec2():
    return io.readspec(data_path('PH957_f.fits'))


@pytest.fixture
def specm(spec,spec2):
    specm = ltsu.collate([spec,spec2])
    return specm


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_readwrite_meta_as_dicts(spec):
    sp = XSpectrum1D.from_tuple((np.array([5,6,7]), np.ones(3), np.ones(3)*0.1))
    sp.meta['headers'][0] = dict(a=1, b='abc')
    sp2 = XSpectrum1D.from_tuple((np.array([8,9,10]), np.ones(3), np.ones(3)*0.1))
    sp2.meta['headers'][0] = dict(c=2, d='efg')
    spec = ltsu.collate([sp,sp2])
    # Write
    spec.write_to_fits(data_path('tmp.fits'))
    spec.write_to_hdf5(data_path('tmp.hdf5'))
    # Read and test
    newspec = io.readspec(data_path('tmp.hdf5'))
    assert newspec.meta['headers'][0]['a'] == 1
    assert newspec.meta['headers'][0]['b'] == 'abc'
    newspec2 = io.readspec(data_path('tmp.fits'))
    assert 'METADATA' in newspec2.meta['headers'][0].keys()

def test_write(spec,specm):
    # FITS
    spec.write(data_path('tmp.fits'))
    spec.write(data_path('tmp.fits'), FITS_TABLE=True)
    # ASCII
    spec.write(data_path('tmp.ascii'))
    # HDF5
    specm.write(data_path('tmp.hdf5'))


def test_hdf5(specm):
    import h5py
    # Write. Should be replaced with tempfile.TemporaryFile
    specm.write_to_hdf5(data_path('tmp.hdf5'))
    #
    specread = io.readspec(data_path('tmp.hdf5'))
    # check a round trip works
    np.testing.assert_allclose(specm.wavelength, specread.wavelength)
    # Add to existing file
    tmp2 = h5py.File(data_path('tmp2.hdf5'), 'w')
    foo = tmp2.create_group('boxcar')
    specm.add_to_hdf5(tmp2, path='/boxcar/')
    tmp2.close()
    # check a round trip works
    spec3 = io.readspec(data_path('tmp2.hdf5'), path='/boxcar/')
    np.testing.assert_allclose(specm.wavelength, spec3.wavelength)


def test_print_repr(spec):
    print(repr(spec))
    print(spec)


def test_write_ascii(spec):
    # Write. Should be replaced with tempfile.TemporaryFile
    spec.write_to_ascii(data_path('tmp.ascii'))
    #
    specb = io.readspec(data_path('tmp.ascii'))
    # check a round trip works
    np.testing.assert_allclose(spec.wavelength, specb.wavelength)


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



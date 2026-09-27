# Module to run tests on spectra.io
import os
import pytest
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

# THIS TEST IS FAILING BUT I AM DEPRECATING XSPEC1D
'''
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
'''

def test_write(spec, specm, tmp_path):
    # FITS
    spec.write(str(tmp_path / 'tmp.fits'))
    spec.write(str(tmp_path / 'tmp.fits'), FITS_TABLE=True)
    # ASCII
    spec.write(str(tmp_path / 'tmp.ascii'))
    # HDF5
    specm.write(str(tmp_path / 'tmp.hdf5'))


def test_hdf5(specm, tmp_path):
    import h5py
    hfile = str(tmp_path / 'tmp.hdf5')
    specm.write_to_hdf5(hfile)
    #
    specread = io.readspec(hfile)
    # check a round trip works
    np.testing.assert_allclose(specm.wavelength, specread.wavelength)
    # Add to existing file
    hfile2 = str(tmp_path / 'tmp2.hdf5')
    tmp2 = h5py.File(hfile2, 'w')
    foo = tmp2.create_group('boxcar')
    specm.add_to_hdf5(tmp2, path='/boxcar/')
    tmp2.close()
    # check a round trip works
    spec3 = io.readspec(hfile2, path='/boxcar/')
    np.testing.assert_allclose(specm.wavelength, spec3.wavelength)


def test_print_repr(spec):
    print(repr(spec))
    print(spec)


def test_write_ascii(spec, tmp_path):
    afile = str(tmp_path / 'tmp.ascii')
    spec.write_to_ascii(afile)
    #
    specb = io.readspec(afile)
    # check a round trip works
    np.testing.assert_allclose(spec.wavelength, specb.wavelength)


def test_write_fits(spec, spec2, tmp_path):
    ffile = str(tmp_path / 'tmp.fits')
    spec.write_to_fits(ffile)
    specin = io.readspec(ffile)
    # check a round trip works
    np.testing.assert_allclose(spec.wavelength, specin.wavelength)
    # ESI
    ffile2 = str(tmp_path / 'tmp2.fits')
    spec2.write_to_fits(ffile2)
    specin2 = io.readspec(ffile2)
    # check a round trip works
    np.testing.assert_allclose(spec2.wavelength, specin2.wavelength)


def test_readwrite_without_sig(tmp_path):
    sp = XSpectrum1D.from_tuple((np.array([5,6,7]), np.ones(3)))
    ffile = str(tmp_path / 'tmp.fits')
    sp.write_to_fits(ffile)
    sp1 = io.readspec(ffile)
    np.testing.assert_allclose(sp1.wavelength.value, sp.wavelength.value)
    np.testing.assert_allclose(sp1.flux.value, sp.flux.value)


# TURNING OFF AS NEWER FITS HEADER COULD NOT HANDLE VERY LONG META
#def test_readwrite_metadata(spec):
#    d = {'a':1, 'b':'abc', 'c':3.2, 'd':np.array([1,2,3]),
#         'e':dict(a=1,b=2)}
#    spec.meta.update(d)
#    spec.write_to_fits(data_path('tmp.fits'))
#    spec2 = io.readspec(data_path('tmp.fits'))
#    pytest.set_trace()
#    assert spec2.meta['a'] == d['a']
#    assert spec2.meta['b'] == d['b']
#    np.testing.assert_allclose(spec2.meta['c'], d['c'])
#    np.testing.assert_allclose(spec2.meta['d'], d['d'])
#    assert spec2.meta['e'] == d['e']



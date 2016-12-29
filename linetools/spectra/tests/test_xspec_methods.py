""" Module to run I/O tests on XSpectrum1D
"""
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import numpy as np
import os
import pytest

import astropy.io.ascii as ascii
from astropy import units as u
import astropy.table

from linetools.spectra import io
from linetools.spectra import utils as lsu
from linetools.spectra.xspectrum1d import XSpectrum1D


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_slice():
    spec = io.readspec(data_path('UM184_nF.fits'))
    spec2 = io.readspec(data_path('PH957_f.fits'))
    spec3 = spec2.copy()
    # Collate to make a multispec spectrum
    mspec = lsu.collate([spec,spec2,spec3])
    # Array
    newspec = mspec[np.array([0,1])]
    # Test
    assert newspec.nspec == 2
    assert not newspec.co_is_set
    # Int
    newspec2 = mspec[1]
    assert newspec2.nspec == 1
    # Slice
    newspec3 = mspec[0:2]
    assert newspec3.nspec == 2


def test_const_sig():
    spec = io.readspec(data_path('UM184_nF.fits'))
    spec.constant_sig(sigv=0.1)
    np.testing.assert_allclose(spec.sig[0], 0.1)


def test_unmask():
    spec = XSpectrum1D.from_file(data_path('UM184_nF.fits'))
    assert np.sum(spec.data['wave'].mask) > 0
    spec.unmask()
    assert np.sum(spec.data['wave'].mask) == 0


def test_addmask():
    spec = XSpectrum1D.from_file(data_path('UM184_nF.fits'))
    assert not spec.data['flux'][0].mask[100]
    mask = spec.data['flux'][0].mask.copy()
    mask[100:110] = True
    spec.add_to_mask(mask)
    assert spec.data['flux'][0].mask[100]
    # Compressed
    badp = spec.flux < 0.1
    spec.add_to_mask(badp, compressed=True)
    assert np.sum(spec.data['flux'][0].mask) > 3000


def test_get_local_s2n():
    spec = XSpectrum1D.from_file(data_path('UM184_nF.fits'))
    wv0 = 4000 * u.AA
    s2n, sig_s2n = spec.get_local_s2n(wv0, 20, flux_th=0.9)
    np.testing.assert_allclose(s2n, 9.30119800567627, rtol=1e-5)
    np.testing.assert_allclose(sig_s2n, 1.0349911451339722, rtol=1e-5)
    # test with continuum
    spec.co = np.ones_like(spec.flux)
    s2n, sig_s2n = spec.get_local_s2n(wv0, 20, flux_th=0.9)
    np.testing.assert_allclose(s2n, 10.330545425415039, rtol=1e-5)
    np.testing.assert_allclose(sig_s2n, 0.4250050187110901, rtol=1e-5)
    # test errors
    # out of range
    with pytest.raises(IOError):
        spec.get_local_s2n(1215*u.AA, 20)
    # sig not defined
    spec = XSpectrum1D.from_tuple((spec.wavelength, spec.flux))
    with pytest.raises(ValueError):
        spec.get_local_s2n(wv0, 20)
    # bad shape for flux_th
    with pytest.raises(ValueError):
        spec.get_local_s2n(wv0, 20, flux_th=np.array([1,2,3,4,5]))
    # npix too big
    with pytest.raises(ValueError):
        spec.get_local_s2n(wv0, 1 + len(spec.wavelength))


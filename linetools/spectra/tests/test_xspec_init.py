# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import numpy as np
import os
import pytest
import pdb

import astropy.io.ascii as ascii
from astropy import units as u
from astropy.io import fits

from specutils.spectrum1d import Spectrum1D

from linetools.spectra import io
from linetools.spectra.xspectrum1d import XSpectrum1D


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

# From arrays
def test_from_tuple():
    idl = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux'],idl['sig']))
    #
    np.testing.assert_allclose(spec.dispersion.value, idl['wave'])
    np.testing.assert_allclose(spec.sig, idl['sig'], atol=2e-3, rtol=0)

    assert spec.dispersion.unit == u.Unit('AA')
    #
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux']))
    np.testing.assert_allclose(spec.dispersion.value, idl['wave'])
    # continuum
    co = np.ones_like(idl['flux'])
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux'],idl['sig'], co))
    np.testing.assert_allclose(spec.dispersion.value, idl['wave'])


# From file
def test_from_spec1d():
    spec = Spectrum1D.from_array(np.array([1,2,3]), np.array([1,1,1]))
    xspec = XSpectrum1D.from_spec1d(spec)

# From file
def test_from_file():
    spec = XSpectrum1D.from_file(data_path('UM184_nF.fits'))
    idl = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])

    np.testing.assert_allclose(spec.dispersion.value, idl['wave'])
    np.testing.assert_allclose(spec.sig, idl['sig'], atol=2e-3, rtol=0)

    assert spec.dispersion.unit == u.Unit('AA')

# Attributes
def test_attrib():
    spec = io.readspec(data_path('UM184_nF.fits'))
    # 
    np.testing.assert_allclose(spec.wvmin.value, 3056.6673905210096)
    np.testing.assert_allclose(spec.wvmax.value, 9205.255609841855)


def test_const_sig():
    spec = io.readspec(data_path('UM184_nF.fits'))
    spec.constant_sig(sigv=0.1)
    np.testing.assert_allclose(spec.sig[0], 0.1)



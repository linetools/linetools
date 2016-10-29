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
from linetools.spectra.xspectrum1d import XSpectrum1D


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_from_tuple():
    tmp = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])
    idl = dict(wave=np.array(tmp['wave']), flux=np.array(tmp['flux']),
               sig=np.array(tmp['sig']))
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux'], idl['sig']))
    #
    np.testing.assert_allclose(spec.data['wave'][spec.select], idl['wave'])
    np.testing.assert_allclose(spec.data['sig'][spec.select], idl['sig'], atol=2e-3, rtol=0)

    assert spec.wavelength.unit == u.Unit('AA')
    #
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux']))
    np.testing.assert_allclose(spec.data['wave'][spec.select], idl['wave'])
    # continuum
    co = np.ones_like(idl['flux'])
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux'],idl['sig'], co))
    np.testing.assert_allclose(spec.data['wave'][spec.select], idl['wave'])

    co = None
    spec = XSpectrum1D.from_tuple((idl['wave'],idl['flux'],idl['sig'], co))
    np.testing.assert_allclose(spec.data['wave'][spec.select], idl['wave'])


def test_from_tuple_Column():
    wv = astropy.table.Column(np.arange(10.), name='wave', unit=None)
    fx = astropy.table.Column(np.ones(len(wv)), name='flux', unit=None)
    sig = np.ones(len(fx))
    spec = XSpectrum1D.from_tuple((wv, fx, sig))
    assert spec.wavelength.unit == u.Angstrom


def test_from_file():
    spec = XSpectrum1D.from_file(data_path('UM184_nF.fits'))
    idl = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])

    np.testing.assert_allclose(spec.data['wave'][spec.select].data, idl['wave'])
    np.testing.assert_allclose(spec.data['sig'][spec.select].data, idl['sig'], atol=2e-3, rtol=0)

    assert spec.wavelength.unit == u.Unit('AA')


def test_attrib():
    spec = io.readspec(data_path('UM184_nF.fits'))
    # 
    np.testing.assert_allclose(spec.wvmin.value, 3056.6673905210096)
    np.testing.assert_allclose(spec.wvmax.value, 9205.255609841855)
    assert spec.npix == 15024


def test_masking():
    wave = 3000. + np.arange(1000)
    flux = np.ones_like(wave)
    sig = 0.1*np.ones_like(wave)
    sig[900:] = 0.
    sig[800:825] = 0.
    wave[900:] = 0.  # WARNING, the data are sorted first!
    #
    spec = XSpectrum1D.from_tuple((wave,flux,sig), masking='edges')
    assert len(spec.wavelength) == 900
    spec2 = XSpectrum1D.from_tuple((wave,flux,sig), masking='all')
    assert len(spec2.wavelength) == 875
    spec3 = XSpectrum1D.from_tuple((wave,flux,sig), masking='none')
    assert len(spec3.wavelength) == len(wave)

def test_co_kludges():
    spec = XSpectrum1D.from_file(data_path('SDSSJ220248.31+123656.3.fits'))
    assert spec.co.size == 4599


def test_errors():

    # from_tuple
    try:
        spec = XSpectrum1D.from_tuple('this_is_not_a_tuple')
    except IOError:
        pass
    try:
        n_tuple = np.array([np.ones(5), np.ones(5)]), np.ones(5)
        spec = XSpectrum1D.from_tuple(n_tuple)
    except IOError:
        pass

    # wrong instances
    flux = [1,2,3]
    wave = [1,2,3]
    try:
        spec = XSpectrum1D(wave, flux)
    except IOError:
        pass

    #wrong shapes
    flux = np.ones(5)
    wave = np.array([1,2,3])
    try:
        spec = XSpectrum1D(wave, flux)
    except IOError:
        pass
    try:
        spec = XSpectrum1D(wave, np.ones(len(wave)), sig=np.ones(2))
    except IOError:
        pass
    try:
        spec = XSpectrum1D(wave, np.ones(len(wave)), co=np.ones(2), verbose = True) # test verbose here too
    except IOError:
        pass

    # wrong masking
    try:
        spec = XSpectrum1D(wave, np.ones(len(wave)), masking = 'wrong_masking')
    except IOError:
        pass

    #wrong units input
    try:
        spec = XSpectrum1D(wave, np.ones(len(wave)), units = 'not_a_dict')
    except IOError:
        pass
    try:
        spec = XSpectrum1D(wave, np.ones(len(wave)), units =dict(wrong_key=2))
    except IOError:
        pass


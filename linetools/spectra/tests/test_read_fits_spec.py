# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals

import os
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.table import Table
import numpy as np
from astropy.io import fits

from linetools.spectra import io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_get_cdelt_dcflag():
    hd = fits.getheader(data_path('UM184_nF.fits'))
    cdelt1, dc_flag = io.get_cdelt_dcflag(hd)
    np.testing.assert_allclose(cdelt1, 3.1870310099e-05)
    np.testing.assert_allclose(dc_flag, 1)


def test_sep_files():
    """ Separate flux/error files"""
    spec = io.readspec(data_path('UM184_nF.fits'))
    idl = ascii.read(data_path('UM184.dat.gz'), names=['wave', 'flux', 'sig'])
    #np.testing.assert_allclose(spec.wavelength.value, idl['wave'])
    np.testing.assert_allclose(spec.data['wave'][spec.select], idl['wave'])
    np.testing.assert_allclose(spec.data['sig'][spec.select], idl['sig'], atol=2e-3, rtol=0)
    assert spec.wavelength.unit == u.Unit('AA')


def test_binary_table():
    spec = io.readspec(data_path('NGC4151sic2a.fits'), head_exten=1)
    # Test head1
    assert len(spec.header['HISTORY']) == 476
    # Data
    np.testing.assert_allclose(spec.flux[0].value, -1.3274571970938327e-16)


def test_setwave():
    hd = fits.getheader(data_path('UM184_nF.fits'))
    wave = io.setwave(hd)
    np.testing.assert_allclose(wave[0], 3042.34515916)
    hd['CRPIX1'] = 10
    wave = io.setwave(hd)
    np.testing.assert_allclose(wave[0], 3040.33648468)


def test_getwave():  # X-Shooter
    spec = io.readspec(data_path('XShooter_XQ100.fits.gz'))
    assert spec.wavelength.unit == u.nm


# ASCII format
def test_read_ascii():
    spec = io.readspec(data_path('UM184.dat.gz'))
    assert spec.wavelength.unit == u.Unit('AA')
    spec = io.readspec(data_path('q0002m422.txt.gz'))
    assert spec.filename.endswith('q0002m422.txt.gz')


def test_uves_popler():
    spec = io.readspec(data_path('popler_sample.fits'))
    assert spec.wavelength.unit == u.Unit('AA')
    assert spec.filename.endswith('popler_sample.fits')


def test_read_table():
    t = Table([(1,2,3), (1,2,3), (1,2,3)], names=['WAVE', 'FLUX', 'ERR'])
    spec = io.readspec(t)
    np.testing.assert_allclose(spec.wavelength[0].value, 1)
    np.testing.assert_allclose(spec.flux[0], 1)
    # Read table with non standard tags
    t = Table([(1,2,3), (1,2,3), (1,2,3)],
              names=['dumb_wave', 'dumb_flux', 'dumb_sig'])
    spec = io.readspec(t, wave_tag='dumb_wave', flux_tag='dumb_flux',
                       sig_tag='dumb_sig')
    np.testing.assert_allclose(spec.wavelength[0].value, 1)
    np.testing.assert_allclose(spec.flux[0], 1)
    # More..
    t = Table([(1,2,3), (1,2,3), (1,2,3)],
              names=['dumb_wave', 'dumb_flux', 'dumb_ivar'])
    spec = io.readspec(t, wave_tag='dumb_wave', flux_tag='dumb_flux',
                       ivar_tag='dumb_ivar')
    t = Table([(1,2,3), (1,2,3), (1,2,3)],
              names=['dumb_wave', 'dumb_flux', 'dumb_var'])
    spec = io.readspec(t, wave_tag='dumb_wave', flux_tag='dumb_flux',
                       var_tag='dumb_var')


def test_errors():
    # no such file
    try:
        io.readspec('filename_that_should_not_exist.txt')
    except IOError:
        pass
    # bad format for ascii with more than 4 columns
    try:
        io.readspec('files/ascii_5columns.txt')
    except IOError:
        pass
    # bad spectra input
    try:
        io.readspec(5) #input as an int
    except IOError:
        pass





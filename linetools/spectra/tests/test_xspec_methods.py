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


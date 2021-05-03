""" Module to run tests on XSpectrum1D attributes
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


@pytest.fixture
def dummy_spec(wave=None, npix=2000, s2n=10., seed=None):
    """
    Parameters
    ----------
    s2n
    seed
    wave

    Returns
    -------
    spec : XSpectrum1D
    """
    if wave is None:
        wave = np.linspace(4000., 5000., npix)
    # Create
    flux = np.ones_like(wave)
    sig = np.ones_like(wave) / s2n
    spec = XSpectrum1D.from_tuple((wave,flux,sig))
    # Noise and append
    if seed is not None:
        spec = spec.add_noise(seed=seed)
    # Return
    return spec


def test_sig():
    dspec = dummy_spec(s2n=10.)
    sig = dspec.sig
    np.testing.assert_allclose(sig[0].value, 0.1)


def test_ivar():
    dspec = dummy_spec(s2n=10.)
    ivar = dspec.ivar
    np.testing.assert_allclose(ivar[0].value, 100)

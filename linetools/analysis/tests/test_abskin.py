from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest
import os
import pdb
from astropy import units as u

from linetools.spectra import io

from linetools.analysis.abskin import generate_stau, pw97_kin


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../../spectra/tests/files')
    return os.path.join(data_dir, filename)


def dummy_spec():
    """ Dummy spectrum for analysis
    """
    wave = np.linspace(1210, 1220, 2000)*u.AA
    npix = wave.size
    fx = np.ones(npix)
    fx[npix//2-5:npix//2+5] = 0.8
    fx[npix//2] = 0.01
    sig = np.ones(npix)*0.1
    wrest = 1215.670*u.AA
    velo = (wave-wrest)/wrest * 3e5*u.km/u.s
    # Return
    return velo, fx, sig


def test_stau():
    # Generate stau
    velo, fx, sig = dummy_spec()
    stau = generate_stau(velo, fx, sig)
    np.testing.assert_allclose(stau[1000], 0.27800134641010432)
    # Different binning now
    stau2 = generate_stau(velo, fx, sig, kbin=35*u.km/u.s)
    np.testing.assert_allclose(stau2[1000], 0.17871515126363854)


def test_pw97():
    # Generate stau
    velo, fx, sig = dummy_spec()
    stau = generate_stau(velo, fx, sig)
    # pw97
    kin_data = pw97_kin(velo, stau)
    np.testing.assert_allclose(kin_data['Dv'].value, 20.)
    np.testing.assert_allclose(kin_data['fedg'], 0.37035142148289424)


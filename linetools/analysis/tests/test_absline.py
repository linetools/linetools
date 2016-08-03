from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest
import pdb
from astropy import units as u

from ..absline import aodm, log_clm, linear_clm, photo_cross, sum_logN

def test_aodm():
    # Make fake spectrum
    wave = np.linspace(1210, 1220, 100)*u.AA
    npix = wave.size
    fx = np.ones(npix)
    fx[npix//2-5:npix//2+5] = 0.8
    fx[npix//2] = 0.01
    sig = np.ones(npix)*0.1
    wrest = 1215.670*u.AA
    fval = 0.4
    velo = (wave-wrest)/wrest * 3e5*u.km/u.s
    # Operate
    N, sig_N, flag_sat = aodm((velo, fx, sig), (wrest, fval))
    # Test
    assert flag_sat is True
    np.testing.assert_allclose((N.value, sig_N.value),
                               (96652191688169.72, 194151305045168.12))
    assert N.unit == u.cm**-2


def test_logclm():
    obj = type(str('Dummy'), (object,), { str('N'): 1e13, str('sig_N'): 5e12 })
    #
    logN, sig_logN = log_clm(obj)
    np.testing.assert_allclose(logN, 13.)


def test_linearclm():
    obj = type(str('Dummy'), (object,), { str('logN'): 13, str('sig_logN'): 0.2 })
    #
    N, sig_N = linear_clm(obj)
    np.testing.assert_allclose(N.value, 1e13)


def test_photocross():
    phto = photo_cross(1, 1, 14.*u.eV)
    assert phto.unit == u.cm**2
    np.testing.assert_allclose(phto.value, 5.870146496955153e-18)


def test_sumlogn_fail():
    obj1 = dict(flag_N=4)
    obj2 = dict(flag_N=4)
    pytest.raises(ValueError, sum_logN, obj1, obj2)
    #
    obj1 = dict(flag_N=1)
    pytest.raises(ValueError, sum_logN, obj1, obj2)


def test_sumlogn_limit():
    obj1 = dict(flag_N=3, logN=15., sig_logN=99.)
    obj2 = dict(flag_N=1, logN=14., sig_logN=0.3)
    flag_N, logN, sig_logN = sum_logN(obj1, obj2)
    # Test
    assert flag_N == 1
    np.testing.assert_allclose((logN, sig_logN), (obj2['logN'], obj2['sig_logN']))


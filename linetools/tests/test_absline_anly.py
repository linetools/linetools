# Module to run tests on AbsLine analysis

# TEST_UNICODE_LITERALS
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio
from linetools import spectralline
from linetools.lists.linelist import LineList


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../spectra/tests/files')
    return os.path.join(data_dir, filename)


def test_aodm_absline():
    # Init CIV 1548
    abslin = AbsLine(1548.195*u.AA, z=2.9304)

    # Set spectrum
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.limits.set([6080.78, 6087.82]*u.AA)
    #abslin.analy['wvlim'] = [6080.78, 6087.82]*u.AA
    #
    abslin.measure_aodm()
    N, sig_N, flgN = [abslin.attrib[key] for key in ['N','sig_N','flag_N']]

    np.testing.assert_allclose(N.value, 76330670518067.16)
    assert N.unit == 1/u.cm**2
    assert flgN == 1

    # Now velocity limits
    abslin.setz(2.92929)
    abslin.limits.set((-150., 150.)*u.km/u.s)
    #
    abslin.measure_aodm()
    N, sig_N, flgN = [abslin.attrib[key] for key in ['N','sig_N','flag_N']]
    np.testing.assert_allclose(N.value, 80369217498895.17)
    return


def test_boxew_absline():
    # Text boxcar EW evaluation
        # Init CIV 1548
    abslin = AbsLine(1548.195*u.AA, z=2.9304)

    # Set spectrum
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.limits.set([6080.78, 6087.82]*u.AA)
    # Measure EW (not rest-frame)
    abslin.measure_ew()
    ew = abslin.attrib['EW']

    np.testing.assert_allclose(ew.value, 0.9935021012055584)
    assert ew.unit == u.AA

    abslin.measure_restew()
    restew = abslin.attrib['EW']
    np.testing.assert_allclose(restew.value, 0.9935021012055584/(1+abslin.z))


def test_gaussew_absline():
    # Text Gaussian EW evaluation
    # Init CIV 1548
    abslin = AbsLine(1548.195*u.AA, z=2.9304)

    # Set spectrum
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.limits.set([6080.78, 6087.82]*u.AA)
    # Measure EW (not rest-frame)
    abslin.measure_ew(flg=2)
    ew = abslin.attrib['EW']

    np.testing.assert_allclose(ew.value, 1.02,atol=0.01)
    assert ew.unit == u.AA

    abslin.measure_ew(flg=2,initial_guesses=(0.5,6081,1))


def test_measurekin_absline():
    # Test Simple kinematics
    abslin = AbsLine('NiII 1741',z=2.307922)

    # Set spectrum
    abslin.analy['spec'] = lsio.readspec(data_path('PH957_f.fits'))
    abslin.limits.set([-70., 70.]*u.km/u.s)

    # Measure Kin
    abslin.measure_kin()
    np.testing.assert_allclose(abslin.attrib['kin']['Dv'].value, 75.)
    np.testing.assert_allclose(abslin.attrib['kin']['fedg'], 0.20005782376000183)


def test_ismatch():
    # Test Simple kinematics
    abslin1 = AbsLine('NiII 1741', z=1.)
    abslin2 = AbsLine('NiII 1741', z=1.)
    # Run
    answer = abslin1.ismatch(abslin2)
    assert answer
    # Tuple too
    answer2 = abslin1.ismatch((1., abslin1.wrest))
    assert answer2


def test_get_tau0():
    abslin1 = AbsLine('HI 1215')
    N = [10**13.0, 10**14.0, 10**20] / (u.cm*u.cm)
    b = [20, 20, 20] * u.km / u.s
    tau0 = abslin1.get_tau0(N, b)


def test_get_Wr_from_N_b():
    abslin1 = AbsLine('HI 1215')
    N = [10**13.0, 10**14.0, 10**20] / (u.cm*u.cm)
    b = [20, 20, 20] * u.km / u.s
    Wr = abslin1.get_Wr_from_N_b(N, b)


def test_get_Wr_from_N_and_viceversa():
    abslin1 = AbsLine('HI 1215')
    N = [10**12.0, 10**12.1, 10**12.2] / (u.cm*u.cm)
    Wr = abslin1.get_Wr_from_N(N)
    N_new = abslin1.get_N_from_Wr(Wr)
    np.testing.assert_allclose(N, N_new, rtol=1e-5)


def test_repr():
    abslin = AbsLine('NiII 1741')
    print(abslin)


def test_manyabslines():
    lines = [1215.670*u.AA, 1025.7222*u.AA, 972.5367*u.AA]*2
    llist = LineList('HI')
    alines = spectralline.many_abslines(lines, llist)


from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools.analysis.voigt import voigt_model 

def test_voigt_sngl_line():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3647))
    # HI line
    abslin = AbsLine(1215.670*u.AA)
    abslin.attrib['N'] = 14.  # log N
    abslin.attrib['b'] = 25.*u.km/u.s
    abslin.attrib['z'] = 2.0
    # Voigt
    vmodel = abslin.generate_voigt(wave=wave)
    np.testing.assert_allclose(vmodel.flux[imn],0.0514565487518)

def test_voigt_multi_line():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3646.2))
    # HI line
    abslin = AbsLine(1215.670*u.AA)
    abslin.attrib['N'] = 17.5  # log N
    abslin.attrib['b'] = 20.*u.km/u.s
    abslin.attrib['z'] = 2.0
    # DI line
    abslin2 = AbsLine('DI 1215')
    abslin2.attrib['N'] = 13.  # log N
    abslin2.attrib['b'] = 15.*u.km/u.s
    abslin2.attrib['z'] = 2.0
    # Voigt
    vmodel3 = voigt_model(wave,[abslin,abslin2])
    np.testing.assert_allclose(vmodel3.flux[imn],0.5714211)

def test_voigt_sngl_tau():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3647))
    # HI line
    abslin = AbsLine(1215.670*u.AA)
    abslin.attrib['N'] = 14.  # log N
    abslin.attrib['b'] = 25.*u.km/u.s
    abslin.attrib['z'] = 2.0
    # Tau
    tau = voigt_model(wave,abslin,flg_ret=2)
    np.testing.assert_allclose(tau[imn], 2.968099279427219)

from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools.analysis.voigt import voigt_from_abslines

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
    np.testing.assert_allclose(vmodel.flux[imn].value,0.05145500775919881)

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
    vmodel3 = voigt_from_abslines(wave,[abslin,abslin2])
    np.testing.assert_allclose(vmodel3.flux[imn].value,0.5715512949324375)

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
    tau = voigt_from_abslines(wave,abslin,ret='tau')
    np.testing.assert_allclose(tau[imn], 2.9681283001576779)

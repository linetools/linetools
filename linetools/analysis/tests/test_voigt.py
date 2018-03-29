from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest
import pdb

from astropy import units as u
from astropy import constants as const

from linetools.spectralline import AbsLine
from linetools.analysis import voigt as lav

c_kms = const.c.to('km/s').value


def test_voigt_model():
    from astropy.modeling import fitting
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    # HI line
    abslin = AbsLine(1215.670*u.AA, z=2.)
    abslin.attrib['N'] = 10**14./u.cm**2
    abslin.attrib['b'] = 25.*u.km/u.s
    # Voigt
    vmodel = abslin.generate_voigt(wave=wave)
    vmodel.sig = 0.1
    # Voigt fit
    abslin.analy['spec'] = vmodel
    abslin.limits.set([-100.,100]*u.km/u.s)
    abslin.measure_aodm(normalize=False)  # Sets analysis pixels
    fitvoigt = lav.single_voigt_model(logN=np.log10(abslin.attrib['N'].value),
        b=abslin.attrib['b'].value, z=2., wrest=abslin.wrest.value,
        gamma=abslin.data['gamma'].value,
        f=abslin.data['f'], fwhm=3.)
    # Restrict parameter space
    fitvoigt.logN.min = 12.
    fitvoigt.b.min = 10.
    fitvoigt.z.min = 2. + -100. * (1 + 2.) / c_kms
    fitvoigt.z.max = 2. + 100 * (1 + 2.) / c_kms

    # Fit
    fitter = fitting.LevMarLSQFitter()
    parm = fitter(fitvoigt,vmodel.wavelength[abslin.analy['pix']].value,
                  vmodel.flux[abslin.analy['pix']].value)
    assert np.abs(parm.logN.value-np.log10(abslin.attrib['N'].value)) < 0.1

def test_voigt_sngl_line():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3647))
    # HI line
    abslin = AbsLine(1215.670*u.AA, z=2.)
    abslin.attrib['N'] = 10**14./u.cm**2
    abslin.attrib['b'] = 25.*u.km/u.s
    # Voigt
    vmodel = abslin.generate_voigt(wave=wave)
    np.testing.assert_allclose(vmodel.flux[imn].value,0.05145500775919881)


def test_voigt_multi_line():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3646.2))
    # HI line
    abslin = AbsLine(1215.670*u.AA, z=2.)
    abslin.attrib['N'] = 10**17.5/u.cm**2
    abslin.attrib['b'] = 20.*u.km/u.s
    # DI line
    abslin2 = AbsLine('DI 1215', z=2.)
    abslin2.attrib['N'] = 10**13./u.cm**2
    abslin2.attrib['b'] = 15.*u.km/u.s
    # Voigt
    vmodel3 = lav.voigt_from_abslines(wave,[abslin,abslin2])
    np.testing.assert_allclose(vmodel3.flux[imn].value,0.5715512949324375)


def test_voigt_fail():
    #
    wave = np.linspace(3644, 3650, 100)
    pytest.raises(ValueError, lav.voigt_from_abslines, wave, [None, None])
    #
    wave = np.linspace(3644, 3650, 100)*u.AA
    pytest.raises(IOError, lav.voigt_from_abslines, wave, [None, None],
                  skip_wveval=True)

def test_voigt_sngl_tau():
    # Wavelength array
    wave = np.linspace(3644, 3650, 100)*u.AA
    imn = np.argmin(np.abs(wave.value-3647))
    # HI line
    abslin = AbsLine(1215.670*u.AA, z=2.)
    abslin.attrib['N'] = 10**14./u.cm**2
    abslin.attrib['b'] = 25.*u.km/u.s
    # Tau
    tau = lav.voigt_from_abslines(wave,abslin,ret='tau')
    assert not np.any(tau<0)
    np.testing.assert_allclose(tau[imn], 2.9681283001576779)


def test_voigt_king():
    vin = np.linspace(0., 1., num=1000)
    a = 0.1
    voigt = lav.voigtking(vin, a)
    # Test
    np.testing.assert_allclose(voigt[50], 0.89440482758173867)


def test_voigt_from_components():
    from linetools.isgm.tests.test_use_abscomp import mk_comp
    wv_array = np.arange(900, 1250, 0.01) * u.AA
    comp1, HIlines = mk_comp('HI', zcomp=0.01, vlim=[-10,10]*u.km/u.s)
    comp2, HIlines = mk_comp('HI', zcomp=0.05, vlim=[-10,10]*u.km/u.s)
    model = lav.voigt_from_components(wv_array, [comp1,comp2])

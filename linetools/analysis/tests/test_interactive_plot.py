import numpy as np
from numpy.random import RandomState

import astropy.units as u
from linetools.analysis import interactive_plot as ltaip
from linetools.spectra.xspectrum1d import XSpectrum1D


def test_local_median():
    fl = np.ones(100)
    wv = np.linspace(1000,2000, 100)
    er = np.ones(100)*0.1
    spec = XSpectrum1D.from_tuple((wv, fl, er))
    lm = ltaip.local_median(spec.wavelength, spec.flux, spec.sig, 1300*u.AA, npix=15)
    np.testing.assert_allclose(lm, 1.)
    # out of ranges
    lm = ltaip.local_median(spec.wavelength, spec.flux, spec.sig, 5300*u.AA, npix=15, default=None)
    assert lm is None
    # bad values
    sig_aux = np.array(spec.sig)
    sig_aux[90:100] = 0.
    spec.sig = sig_aux
    lm = ltaip.local_median(spec.wavelength, spec.flux, spec.sig, 1950*u.AA, npix=2, default=None)
    assert lm is None
    # with real noise now
    rstate = RandomState(3)
    spec2 = spec.add_noise(s2n=spec.flux/spec.sig, rstate=rstate)
    lm = ltaip.local_median(spec2.wavelength, spec2.flux, spec2.sig, 1300*u.AA, npix=15)
    np.testing.assert_allclose(lm, 0.9428995251655579)



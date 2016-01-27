from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np
from astropy.io import fits

from linetools.spectra.plotting import get_flux_plotrange

def test_get_flux_plotrange():
    flux = [1, 2, 0.5, np.nan, 0.7, 1.1]
    qflux = u.quantity.Quantity(flux)
    ref = -0.273, 2.73
    np.testing.assert_allclose(get_flux_plotrange(flux), ref)
    np.testing.assert_allclose(get_flux_plotrange(qflux), ref)

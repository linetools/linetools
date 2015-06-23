# Module to run tests on AbsLine analysis

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../spectra/tests/files')
    return os.path.join(data_dir, filename)

def test_boxew_absline():
	# Init CIV 1548
    abslin = AbsLine(1548.195*u.AA)

    # Set spectrum 
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.analy['WVMNX'] = [6080.78, 6087.82]*u.AA
    ew, sigew = abslin.box_ew() 
    np.testing.assert_allclose(ew, 0.990466303607004*u.AA)

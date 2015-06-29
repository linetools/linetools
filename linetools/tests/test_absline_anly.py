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

'''
def test_aodm_absline():
    # Init CIV 1548
    abslin = AbsLine(1548.195*u.AA)

    # Set spectrum 
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.analy['wvlim'] = [6080.78, 6087.82]*u.AA
    # 
    flg_sat = None
    N, sigN = abslin.aodm(flg_sat=flg_sat) 

    np.testing.assert_allclose(N.value, 299143539014510.06)
    assert N.unit == 1/u.cm**2
'''

def test_boxew_absline():
	# Init CIV 1548
    abslin = AbsLine(1548.195*u.AA)

    # Set spectrum 
    abslin.analy['spec'] = lsio.readspec(data_path('UM184_nF.fits')) # Fumagalli+13 MagE spectrum
    abslin.analy['wvlim'] = [6080.78, 6087.82]*u.AA
    ew, sigew = abslin.box_ew() 

    np.testing.assert_allclose(ew.value, 0.990466303607004)
    assert ew.unit == u.AA

def test_ismatch():
    # Init CIV 1548
    abslin = AbsLine(1548.195*u.AA)
    abslin.attrib['z'] = 1.2322

    abslin2 = AbsLine(1548.195*u.AA)
    abslin2.attrib['z'] = 1.2322

    assert abslin.ismatch(abslin2)

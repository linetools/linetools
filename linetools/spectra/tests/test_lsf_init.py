# Module to run tests on spectra.lsf
import os
import pytest
from astropy import units as u
import numpy as np

from linetools.spectra.lsf import LSF

def test_lsf_COS():
    wv_array = np.arange(1200,1500,1)*u.AA
    
    gratings = ['G130M','G160M', 'G140L','G230L', 'G185M', 'G225M', 'G285M']
    life_positions = ['1']
    
    for grating in gratings:
        for lp in life_positions:
            print(grating,lp)
            instr_config = dict(name='COS',grating=grating,life_position=lp)
            lsf = LSF(wv_array,instr_config)


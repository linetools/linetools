# Module to run tests on spectra.lsf
import os
import pytest
from astropy import units as u
import numpy as np

from linetools.spectra.lsf import LSF

def test_lsf_COS():
    
    gratings = ['G130M','G160M', 'G140L','G230L', 'G185M', 'G225M', 'G285M']
    life_positions = ['1','2']
    cen_waves_G160M = ['1577','1589','1600','1611','1623']
    cen_waves_G130M = ['1291','1300','1309','1318','1327']

    for grating in gratings:
        for lp in life_positions:
            
            instr_config = dict(name='COS',grating=grating,life_position=lp)
            if lp == '2':
            	if grating not in ['G130M','G160M']:
            		continue
            	if grating == 'G130M':
            		cen_waves_aux = cen_waves_G130M
            	elif grating == 'G160M':
            		cen_waves_aux = cen_waves_G160M

            	for cen_wave in cen_waves_aux:
            		instr_config['cen_wave'] = cen_wave
            		lsf = LSF(instr_config)
            		print(lp,grating,cen_wave)
            elif lp == '1':
				lsf = LSF(instr_config)
				print(lp,grating)
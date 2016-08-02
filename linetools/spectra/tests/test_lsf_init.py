# Module to run tests on spectra.lsf

from __future__ import print_function, absolute_import, \
     division, unicode_literals

import os
import pytest
from astropy import units as u
import numpy as np

from linetools.spectra.lsf import LSF

def test_lsf_COS():
    
    gratings = ['G130M','G160M', 'G140L','G230L', 'G185M', 'G225M', 'G285M']
    life_positions = ['1','2','3']
    cen_waves_G160M = ['1577','1589','1600','1611','1623']
    cen_waves_G130M = ['1291','1300','1309','1318','1327']

    for grating in gratings:
        for lp in life_positions:
            
            instr_config = dict(name='COS',grating=grating,life_position=lp)
            if lp in ['2','3']:
                if grating not in ['G130M','G160M']:
                        continue
                if grating == 'G130M':
                        cen_waves_aux = cen_waves_G130M
                elif grating == 'G160M':
                        cen_waves_aux = cen_waves_G160M

                for cen_wave in cen_waves_aux:
                        instr_config['cen_wave'] = cen_wave
                        lsf = LSF(instr_config)
                        print(lp, grating, cen_wave)
            elif lp == '1':
                lsf = LSF(instr_config)
                print(lp, grating)


def test_lsf_init_errors():
    with pytest.raises(TypeError):
        lsf = LSF('not_a_dict')
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(wrong_key='xx'))
    with pytest.raises(NotImplementedError):
        lsf = LSF(dict(name='not_COS'))
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(name='COS', not_grating_given='xx'))
    with pytest.raises(NotImplementedError):
        lsf = LSF(dict(name='COS', grating='not_implemented_grating'))
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(name='COS', grating='G130M', not_life_pos_given='xx'))
    with pytest.raises(ValueError):
        lsf = LSF(dict(name='COS', grating='G130M', life_position='-1'))
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(name='COS', grating='G130M', life_position='2', no_cen_wave_given='xx'))


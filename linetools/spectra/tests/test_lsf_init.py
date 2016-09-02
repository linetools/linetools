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


def test_lsf_STIS():
    gratings = ['G750L', 'G140M', 'G140L', 'G230M', 'G230L']
                # Echelle gratings, G750M, G430M, G430L, G230LB G230MB still not fully implemented.
                #todo : implement Echelle gratings, G430M, G430L, G230LB G230MB

    available_slits = {
            'G140L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G140M': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G230L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G230M': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'E140H': ['0.1x0.03', '0.2x0.09', '0.2x0.2', '6x0.2'],
            'E140M': ['0.1x0.03', '0.2x0.06', '0.2x0.2', '6x0.2'],
            'E230H': ['0.1x0.03', '0.1x0.09', '0.1x0.2', '6x0.2'],
            'E230M': ['0.1x0.03', '0.2x0.06', '0.2x0.2', '6x0.2'],
            'G430L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G750L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0']
    }

    for grating in gratings:
        # slits
        for slit in available_slits[grating]:
            instr_config = dict(name='STIS',grating=grating, slit=slit)
            lsf = LSF(instr_config)


def test_lsf_init_errors():
    with pytest.raises(TypeError):
        lsf = LSF('not_a_dict')
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(wrong_key='xx'))
    with pytest.raises(NotImplementedError):
        lsf = LSF(dict(name='bad_instrument'))

    # for COS
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

    # for STIS
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(name='STIS', not_grating_given='xx'))
    with pytest.raises(NotImplementedError):
        lsf = LSF(dict(name='STIS', grating='not_implemented_grating'))
    with pytest.raises(SyntaxError):
        lsf = LSF(dict(name='STIS', grating='G140L', not_slit_given='xx'))
    with pytest.raises(NotImplementedError):
        lsf = LSF(dict(name='STIS', grating='G140L', slit='bad_slit'))

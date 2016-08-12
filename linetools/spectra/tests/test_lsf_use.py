# Module to run tests on spectra.lsf
from __future__ import print_function, absolute_import, \
     division, unicode_literals


import os
import pytest
import astropy.units as u
import numpy as np
from astropy.table import Table

from linetools.spectra.lsf import LSF

def test_interpolate_to_wv0(plot=False):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv0()'
    wv0 = 1160*u.AA
    cos_dict = dict(name='COS',grating='G130M',life_position='2',cen_wave='1309')
    lsf_cos = LSF(cos_dict)
    lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
    assert lsf_tab[len(lsf_tab)//2]['wv'] == wv0.value, err_msg
    assert lsf_tab[len(lsf_tab)//2]['kernel'] == np.max(lsf_tab['kernel']), err_msg
    if plot:
        import matplotlib.pyplot as plt
        wv_array = np.arange(1200,1400,10)*u.AA
        for wv in wv_array:
            lsf_tab = lsf_cos.interpolate_to_wv0(wv)
            plt.plot(lsf_tab['wv']-wv.value,lsf_tab['kernel'],'-')
        plt.show()

def test_interpolate_to_wv0_wv0shortlong(plot=False):
    err_msg = 'Something is wrong with short wavelength handling in LSF.interpolate_to_wv0()'
    wv0 = 1105.0 * u.AA
    cos_dict = dict(name='COS', grating='G130M', life_position='1')
    lsf_cos = LSF(cos_dict)
    lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
    assert lsf_tab[len(lsf_tab) // 2]['wv'] == wv0.value, err_msg
    assert lsf_tab[len(lsf_tab) // 2]['kernel'] == np.max(lsf_tab['kernel']), err_msg
    if plot:
        import matplotlib.pyplot as plt
        lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
        plt.plot(lsf_tab['wv'] - wv0.value, lsf_tab['kernel'], '-')
        plt.show()
    err_msg = 'Something is wrong with long wavelength handling in LSF.interpolate_to_wv0()'
    wv0 = 1796.0 * u.AA
    cos_dict = dict(name='COS', grating='G160M', life_position='1')
    lsf_cos = LSF(cos_dict)
    lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
    assert lsf_tab[len(lsf_tab) // 2]['wv'] == wv0.value, err_msg
    assert lsf_tab[len(lsf_tab) // 2]['kernel'] == np.max(lsf_tab['kernel']), err_msg
    if plot:
        import matplotlib.pyplot as plt
        lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
        plt.plot(lsf_tab['wv'] - wv0.value, lsf_tab['kernel'], '-')
        plt.show()

def test_interpolate_to_wv_array(plot=False):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv_array()'
    wv_array = np.arange(1600,1601,0.001)*u.AA
    wv_array = np.arange(1600,1650,0.001)*u.AA
    cen_waves = ['1577','1589A','1600','1611','1623']
    colors = ['k','b','g','r','orange']
    lsf_dict = dict()
    for i,cen_wave in enumerate(cen_waves):
        cos_dict_aux = dict(name='COS',grating='G160M',life_position='2',cen_wave=cen_wave)
        lsf_dict[cen_wave] = LSF(cos_dict_aux)
        lsf_tab = lsf_dict[cen_wave].interpolate_to_wv_array(wv_array)
        assert isinstance(lsf_tab,Table), err_msg
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wv_array,lsf_tab['kernel'],'-',color=colors[i])
    if plot:
        plt.show()

def test_get_lsf(plot=False):
    err_msg = 'Something is wrong with LSF.get_lsf()'
    wv_array = np.arange(1250,1251,0.0001)*u.AA
    cen_waves = ['1291','1300','1309','1318A','1327A']
    colors = ['k','b','g','r','orange']
    lsf_dict = dict()
    for i,cen_wave in enumerate(cen_waves):
        cos_dict_aux = dict(name='COS',grating='G130M',life_position='2',cen_wave=cen_wave)
        lsf_dict[cen_wave] = LSF(cos_dict_aux)
        lsf_kernel = lsf_dict[cen_wave].get_lsf(wv_array)
        assert isinstance(lsf_kernel,np.ndarray), err_msg
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wv_array,lsf_kernel,'-',color=colors[i])
    if plot:
        plt.show()



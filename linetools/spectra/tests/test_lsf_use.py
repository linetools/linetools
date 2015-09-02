# Module to run tests on spectra.lsf
import os
import pytest
from astropy import units as u
import numpy as np
# import matplotlib.pyplot as plt

from linetools.spectra.lsf import LSF

def test_interpolate_to_wv0(plot=False):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv0()'
    wv0 = 1160*u.AA
    cos_dict = dict(name='COS',grating='G130M',life_position='2',cen_wave='1309')
    lsf_cos = LSF(cos_dict)
    lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
    assert lsf_tab[len(lsf_tab)/2]['wv'] == wv0.value, err_msg
    assert lsf_tab[len(lsf_tab)/2]['kernel'] == np.max(lsf_tab['kernel']), err_msg
    if plot:
        wv_array = np.arange(1200,1400,10)*u.AA
        for wv in wv_array:
            lsf_tab = lsf_cos.interpolate_to_wv0(wv)
            plt.plot(lsf_tab['wv']-wv.value,lsf_tab['kernel'],'-')
        plt.show()

def test_interpolate_to_wv_array(plot=False):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv_array()'
    wv_array = np.arange(1600,1601,0.001)*u.AA
    cen_waves = ['1577','1589','1600','1611','1623']
    colors = ['k','b','g','r','orange']
    lsf_dict = dict()
    for i,cen_wave in enumerate(cen_waves):
        cos_dict_aux = dict(name='COS',grating='G160M',life_position='2',cen_wave=cen_wave)
        lsf_dict[cen_wave] = LSF(cos_dict_aux)
        if plot:
            plt.plot(wv_array,lsf_dict[cen_wave].interpolate_to_wv_array(wv_array),'-',color=colors[i])
    if plot:
        plt.show()



# Module to run tests on spectra.lsf
from __future__ import print_function, absolute_import, \
    division, unicode_literals

import os
import pytest
import astropy.units as u
import numpy as np
from astropy.table import Table

from linetools.spectra.lsf import LSF


def test_interpolate_to_wv0(plot=False, lp='2'):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv0()'
    wv0 = 1160 * u.AA
    cos_dict = dict(name='COS', grating='G130M', life_position=lp, cen_wave='1309')
    stis_dict = dict(name='STIS', grating='G140L', slit='52x0.2')
    lsf_cos = LSF(cos_dict)
    lsf_stis = LSF(stis_dict)
    for lsf in [lsf_cos, lsf_stis]:
        lsf_tab = lsf.interpolate_to_wv0(wv0)
        assert lsf_tab[len(lsf_tab) // 2]['wv'] == wv0.value, err_msg
        assert lsf_tab[len(lsf_tab) // 2]['kernel'] == np.max(lsf_tab['kernel']), err_msg
        if plot:
            import matplotlib.pyplot as plt
            if lp == '3':
                wv_array = np.linspace(1200, 1400, 10) * u.AA
            else:
                wv_array = np.arange(1200, 1400, 10) * u.AA

            for wv in wv_array:
                lsf_tab = lsf_cos.interpolate_to_wv0(wv)
                plt.plot(lsf_tab['wv'] - wv.value, lsf_tab['kernel'], '-')
                plt.suptitle(lsf.name)
            # import pdb; pdb.set_trace()
            plt.show()
    # test last column interpolation
    lsf = LSF(dict(name='COS', grating='G130M', life_position='1'))
    tab = lsf.interpolate_to_wv0(1450 * u.AA)


def test_interpolate_to_wv_array(plot=False, lp='2'):
    err_msg = 'Something is wrong with LSF.interpolate_to_wv_array()'
    wv_array = np.arange(1600, 1601, 0.001) * u.AA
    wv_array = np.arange(1600, 1650, 0.001) * u.AA
    cen_waves = ['1577', '1589A', '1600', '1611', '1623']
    colors = ['k', 'b', 'g', 'r', 'orange']
    lsf_dict = dict()
    for i, cen_wave in enumerate(cen_waves):
        cos_dict_aux = dict(name='COS', grating='G160M', life_position=lp, cen_wave=cen_wave)
        lsf_dict[cen_wave] = LSF(cos_dict_aux)
        lsf_tab = lsf_dict[cen_wave].interpolate_to_wv_array(wv_array)
        assert isinstance(lsf_tab, Table), err_msg
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wv_array, lsf_tab['kernel'], '-', color=colors[i])
    if plot:
        plt.show()
    # other tests
    lsf = LSF(dict(name='COS', grating='G130M', life_position='1'))
    wv_array = np.linspace(1210, 1211, 100) * u.AA
    # cubic
    tab = lsf.interpolate_to_wv_array(wv_array, kind='cubic', debug=True)

    # errors
    with pytest.raises(SyntaxError):
        tbl = lsf.interpolate_to_wv_array('bad_input')
    with pytest.raises(SyntaxError):
        x = np.array([[1, 2, 3], [4, 5, 6]])
        tbl = lsf.interpolate_to_wv_array(x)  # bad shape
    with pytest.raises(ValueError):
        tbl = lsf.interpolate_to_wv_array(np.array([1, 2]) * u.AA, kind='wrong_kind')
    with pytest.raises(ValueError):
        tbl = lsf.interpolate_to_wv_array(np.array([1, 2]) * u.AA, kind='cubic')  # bad input wv_array

def test_get_lsf(plot=False, lp='2'):
    err_msg = 'Something is wrong with LSF.get_lsf()'
    wv_array = np.arange(1250, 1251, 0.0001) * u.AA
    cen_waves = ['1291', '1300', '1309', '1318A', '1327A']
    colors = ['k', 'b', 'g', 'r', 'orange']
    lsf_dict = dict()
    for i, cen_wave in enumerate(cen_waves):
        cos_dict_aux = dict(name='COS', grating='G130M', life_position=lp, cen_wave=cen_wave)
        lsf_dict[cen_wave] = LSF(cos_dict_aux)
        lsf_kernel = lsf_dict[cen_wave].get_lsf(wv_array)
        assert isinstance(lsf_kernel, np.ndarray), err_msg
        if plot:
            import matplotlib.pyplot as plt
            plt.plot(wv_array, lsf_kernel, '-', color=colors[i])
    if plot:
        plt.show()


def test_all_lp(plot=False):
    for lp in ['1', '3']:  # lp='2' will be tested as default
        test_interpolate_to_wv0(plot=plot, lp=lp)
        test_interpolate_to_wv_array(plot=plot, lp=lp)
        test_get_lsf(plot=plot, lp=lp)


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
    # test error outside range
    wv0 = 1300.0 * u.AA
    with pytest.raises(ValueError):
        lsf_tab = lsf_cos.interpolate_to_wv0(wv0)
    wv0 = 1900.0 * u.AA
    with pytest.raises(ValueError):
        lsf_tab = lsf_cos.interpolate_to_wv0(wv0)

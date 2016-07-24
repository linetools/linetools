# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import numpy as np
import pdb

#import matplotlib
#matplotlib.use('Agg')

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio
from linetools.analysis import absline as ltaa
from linetools.isgm import utils as ltiu
import linetools.utils as ltu

import imp
lt_path = imp.find_module('linetools')[1]

#import pdb
#pdb.set_trace()
# Set of Input lines

def mk_comp(ctype,vlim=[-300.,300]*u.km/u.s,add_spec=False, use_rand=True,
            add_trans=False, zcomp=2.92939):
    # Read a spectrum Spec
    if add_spec:
        xspec = lsio.readspec(lt_path+'/spectra/tests/files/UM184_nF.fits')
    else:
        xspec = None
    # AbsLines
    if ctype == 'HI':
        all_trans = ['HI 1215', 'HI 1025']
    elif ctype == 'SiII':
        all_trans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
        if add_trans:
            all_trans += ['SiII 1193']
    abslines = []
    for trans in all_trans:
        iline = AbsLine(trans)
        iline.attrib['z'] = zcomp
        if use_rand:
            rnd = np.random.rand()
        else:
            rnd = 0.
        iline.attrib['logN'] = 13.3 + rnd
        iline.attrib['sig_logN'] = 0.15
        iline.attrib['flag_N'] = 1
        iline.analy['spec'] = xspec
        iline.analy['vlim'] = vlim
        iline.analy['wvlim'] = iline.wrest * (1 + zcomp + ltu.give_dz(vlim, zcomp))
        _,_ = ltaa.linear_clm(iline.attrib)  # Loads N, sig_N
        abslines.append(iline)
    # Component
    abscomp = AbsComponent.from_abslines(abslines)
    return abscomp, abslines


def test_add_absline():
    abscomp,_ = mk_comp('HI', zcomp=0)
    abscomp.add_absline(AbsLine('HI 972'), chk_sep=False, chk_vel=False)
    with pytest.raises(ValueError):
        abscomp.add_absline(AbsLine('HI 949'), vtoler=-10)
    # failed addition
    bad_absline = AbsLine('CIV 1550')
    bad_absline.analy['vlim'] = [500, 1000]*u.km/u.s
    bad_absline.attrib['coord'] = SkyCoord(20,20, unit='deg')
    abscomp.add_absline(bad_absline)


def test_fromtodict():
    SiIIcomp1,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    cdict = SiIIcomp1.to_dict()
    #
    assert isinstance(cdict, dict)
    assert cdict['Zion'] == (14, 2)
    # And instantiate
    newcomp = AbsComponent.from_dict(cdict)
    assert isinstance(newcomp, AbsComponent)
    newcomp = AbsComponent.from_dict(cdict, coord=SkyCoord(0,0, unit='deg'))


def test_build_table():
    abscomp,_ = mk_comp('HI')
    # Instantiate
    comp_tbl = abscomp.build_table()
    # Test
    assert isinstance(comp_tbl,QTable)
    # empty
    abscomp._abslines = []
    comp_tbl = abscomp.build_table()


def test_synthesize_colm():
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True,
                        add_trans=True)
    # Column
    abscomp.synthesize_colm(redo_aodm=True)
    # Test
    np.testing.assert_allclose(abscomp.logN, 13.594445560856554)
    # Reset flags (for testing)
    abscomp2,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True, use_rand=False,
                         add_trans=True)
    for iline in abscomp2._abslines:
        if iline.data['name'] == 'SiII 1260':
            iline.attrib['flag_N'] = 2
        elif iline.data['name'] == 'SiII 1808':
            iline.attrib['flag_N'] = 3
        elif iline.data['name'] == 'SiII 1193':
            iline.attrib['flag_N'] = 0
        else:
            iline.attrib['flag_N'] = 1
    abscomp2.synthesize_colm()
    # Test
    np.testing.assert_allclose(abscomp2.logN, 13.3)
    # Another
    abscomp3,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True, use_rand=False)
    for iline in abscomp3._abslines:
        if iline.data['name'] == 'SiII 1260':
            iline.attrib['flag_N'] = 3
        elif iline.data['name'] == 'SiII 1808':
            iline.attrib['flag_N'] = 2
        else:
            iline.attrib['flag_N'] = 1
    abscomp3.synthesize_colm()
    # Test
    np.testing.assert_allclose(abscomp3.logN, 13.3)
    # test error
    with pytest.raises(IOError):
        abscomp3.synthesize_colm(overwrite=False)
    with pytest.raises(ValueError):
        abscomp3._abslines[0].attrib['N'] = 0 / u.cm / u.cm
        abscomp3.synthesize_colm(overwrite=True)

def test_build_components_from_lines():
    # Lines
    abscomp,HIlines = mk_comp('HI')
    abscomp,SiIIlines = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # Components
    comps = ltiu.build_components_from_abslines([HIlines[0],HIlines[1],SiIIlines[0],SiIIlines[1]])
    assert len(comps) == 2


def test_iontable_from_components():
    # Lines
    abscomp,HIlines = mk_comp('HI')
    abscomp,SiIIlines = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # Components
    comps = ltiu.build_components_from_abslines([HIlines[0],HIlines[1],SiIIlines[0],SiIIlines[1]])
    tbl = ltiu.iontable_from_components(comps)
    assert len(tbl) == 2


def test_cog():
    # Component
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # COG
    COG_dict = abscomp.cog(redo_EW=True)
    # Test
    np.testing.assert_allclose(COG_dict['logN'],13.693355878125537)
    np.testing.assert_allclose(COG_dict['sig_logN'],0.054323725737309987)


def test_synthesize_components():
    #
    SiIIcomp1,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    SiIIcomp1.synthesize_colm(redo_aodm=True)
    #
    SiIIcomp2,_ = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, add_spec=True)
    SiIIcomp2.synthesize_colm(redo_aodm=True)
    #
    synth_SiII = ltiu.synthesize_components([SiIIcomp1,SiIIcomp2])
    np.testing.assert_allclose(synth_SiII.logN, 13.862454764546792)
    np.testing.assert_allclose(synth_SiII.sig_logN, 0.010146946475971825)
    # Failures
    pytest.raises(IOError, ltiu.synthesize_components, 1)
    pytest.raises(IOError, ltiu.synthesize_components, [1,2])

"""
def test_stack_plot():
    # AbsLine(s)
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    # Spectra
    xspec = lsio.readspec(lt_path+'/spectra/tests/files/UM184_nF.fits')
    lya.analy['spec'] = xspec
    lyb.analy['spec'] = xspec
    # Instantiate
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # Plot
    abscomp.stack_plot(show=False)
"""


def test_repr_vpfit():
    abscomp, HIlines = mk_comp('HI')
    s = abscomp.repr_vpfit()
    assert s == 'HI 2.92939 0.00000 10.00 0.00 0.00 0.00\n'

    s = abscomp.repr_vpfit(b=15*u.km/u.s)
    assert s == 'HI 2.92939 0.00000 15.00 0.00 0.00 0.00\n'

    abscomp.comment = 'Something'
    s = abscomp.repr_vpfit()
    assert s == 'HI 2.92939 0.00000 10.00 0.00 0.00 0.00! Something\n'
    s = abscomp.repr_vpfit(tie_strs=('a', 'b', 'CD'), fix_strs=('', 'f', ''))
    assert s == 'HI 2.92939a 0.00000 10.00F 0.00 0.00cd 0.00! Something\n'

    abscomp, SiIIlines = mk_comp('SiII')
    s = abscomp.repr_vpfit()
    assert s == 'SiII 2.92939 0.00000 10.00 0.00 0.00 0.00\n'

    # errors
    with pytest.raises(TypeError):
        s = abscomp.repr_vpfit(tie_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_vpfit(fix_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_vpfit(fix_strs=('1','2','3','4','5'))



def test_repr_alis():
    abscomp, HIlines = mk_comp('HI')
    s = abscomp.repr_alis()
    assert s == 'voigt   ion=1H_I 0.00 redshift=2.92939 0.0 1.0E+04\n'

    abscomp, SiIIlines = mk_comp('SiII')
    s = abscomp.repr_alis(T_kin=10**5*u.K)
    assert s == 'voigt   ion=28Si_II 0.00 redshift=2.92939 0.0 1.0E+05\n'

    abscomp.comment = 'Something'
    s = abscomp.repr_alis()
    assert s == 'voigt   ion=28Si_II 0.00 redshift=2.92939 0.0 1.0E+04# Something\n'
    s = abscomp.repr_alis(tie_strs=('a', 'b', 'CD',''), fix_strs=('', 'f', '', ''))
    assert s == 'voigt   ion=28Si_II 0.00a redshift=2.92939F 0.0cd 1.0E+04# Something\n'

    # errors
    with pytest.raises(TypeError):
        s = abscomp.repr_alis(tie_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_alis(fix_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_alis(fix_strs=('1','2','3','4','5'))


def test_get_wvobs_chunks():
    abscomp, HIlines = mk_comp('HI', zcomp=0, vlim=[0,10]*u.km/u.s)
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    np.testing.assert_allclose(wvobs_chunks[0][0], 1215.67*u.AA)
    np.testing.assert_allclose(wvobs_chunks[0][1], 1215.71055106*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][0], 1025.7222*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][1], 1025.75641498*u.AA)
    abscomp, HIlines = mk_comp('HI', zcomp=1, vlim=[-100,100]*u.km/u.s)
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    np.testing.assert_allclose(wvobs_chunks[0][0], 2430.52912749*u.AA)
    np.testing.assert_allclose(wvobs_chunks[0][1], 2432.15114303*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][0], 2050.76022589*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][1], 2052.12880236*u.AA)


def test_coincident_components():
    abscomp, HIlines = mk_comp('HI', zcomp=2.92939)
    SiIIcomp1,_ = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, zcomp=2.92939)
    SiIIcomp2,_ = mk_comp('SiII',vlim=[-300.,0.]*u.km/u.s, zcomp=2.92939)
    assert ltiu.coincident_components(abscomp, abscomp)  # should overlap
    assert not ltiu.coincident_components(abscomp, SiIIcomp1)  # should not overlap
    assert not ltiu.coincident_components(SiIIcomp2, SiIIcomp1) # should not overlap

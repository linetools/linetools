# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.table import QTable
import numpy as np
import pdb

#import matplotlib
#matplotlib.use('Agg')

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio
from linetools.analysis import absline as ltaa
from linetools.isgm import utils as ltiu

import imp
lt_path = imp.find_module('linetools')[1]

#import pdb
#pdb.set_trace()
# Set of Input lines

def mk_comp(ctype,vlim=[-300.,300]*u.km/u.s,add_spec=False):
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
    abslines = []
    for trans in all_trans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.attrib['logN'] = 13.3 + np.random.rand()
        iline.attrib['sig_logN'] = 0.15
        iline.attrib['flag_N'] = 1
        iline.analy['spec'] = xspec
        iline.analy['vlim'] = vlim
        _,_ = ltaa.linear_clm(iline.attrib) # Loads N, sig_N
        abslines.append(iline)
    # Component
    abscomp = AbsComponent.from_abslines(abslines)
    return abscomp

def test_build_table():
    abscomp = mk_comp('HI')
    # Instantiate
    comp_tbl = abscomp.build_table()
    # Test
    assert isinstance(comp_tbl,QTable)

def test_synthesize_colm():
    abscomp = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # Column
    abscomp.synthesize_colm(redo_aodm=True)
    # Test
    np.testing.assert_allclose(abscomp.logN,13.594447075294818)

def test_cog():
    # Component
    abscomp = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # COG
    COG_dict = abscomp.cog(redo_EW=True)
    # Test
    np.testing.assert_allclose(COG_dict['logN'],13.693355878125537)
    np.testing.assert_allclose(COG_dict['sig_logN'],0.054323725737309987)

def test_synthesize_components():
    #
    SiIIcomp1 = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    SiIIcomp1.synthesize_colm(redo_aodm=True)
    #
    SiIIcomp2 = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, add_spec=True)
    SiIIcomp2.synthesize_colm(redo_aodm=True)
    #
    synth_SiII = ltiu.synthesize_components([SiIIcomp1,SiIIcomp2])
    np.testing.assert_allclose(synth_SiII.logN,13.862456155250918)
    np.testing.assert_allclose(synth_SiII.sig_logN,0.010146948602759272)
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

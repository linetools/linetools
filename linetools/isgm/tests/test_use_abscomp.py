# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.table import QTable
import numpy as np

#import matplotlib
#matplotlib.use('Agg')

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio

import imp
lt_path = imp.find_module('linetools')[1]

#import pdb
#pdb.set_trace()
# Set of Input lines

def test_build_table():
    # AbsLine(s)
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    # Instantiate
    abscomp = AbsComponent.from_abslines([lya,lyb])
    comp_tbl = abscomp.build_table()
    # Test
    assert isinstance(comp_tbl,QTable)


def test_synthesize_colm():
    # Read a spectrum Spec
    xspec = lsio.readspec(lt_path+'/spectra/tests/files/UM184_nF.fits')
    # AbsLines
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        iline.analy['spec'] = xspec
        abslines.append(iline)
    # Component
    abscomp = AbsComponent.from_abslines(abslines)
    # Column
    abscomp.synthesize_colm(redo_aodm=True)
    # Test
    np.testing.assert_allclose(abscomp.logN,13.594447075294818)

def test_cog():
    # Read a spectrum Spec
    xspec = lsio.readspec(lt_path+'/spectra/tests/files/UM184_nF.fits')
    # AbsLines
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        iline.analy['spec'] = xspec
        abslines.append(iline)
    # Component
    abscomp = AbsComponent.from_abslines(abslines)
    # COG
    COG_dict = abscomp.cog(redo_EW=True)
    # Test
    np.testing.assert_allclose(COG_dict['logN'],13.693355878125537)
    np.testing.assert_allclose(COG_dict['sig_logN'],0.054323725737309987)

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

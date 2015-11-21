# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.table import QTable
import numpy as np

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.spectra import io as lsio

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

def test_stack_plot():
    import imp
    lt_path = imp.find_module('linetools')[1]
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

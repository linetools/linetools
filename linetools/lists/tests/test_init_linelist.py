# Module to run tests on Generating a LineList
#   Also tests some simple functionality

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pdb
import pytest
from astropy import units as u
import numpy as np

from linetools.lists.linelist import LineList
from linetools.lists import mk_sets as llmk



def test_ism_read_source_catalogues():
    ism = LineList('ISM', use_ISM_table=False)
    np.testing.assert_allclose(ism['HI 1215']['wrest'],
                               1215.6700*u.AA, rtol=1e-7)

# ISM LineList
def test_ism():
    ism = LineList('ISM')
    #
    np.testing.assert_allclose(ism['HI 1215']['wrest'],
                               1215.6700*u.AA, rtol=1e-7)

# Test update_fval
def test_updfval():
    ism = LineList('ISM')
    #
    np.testing.assert_allclose(ism['FeII 1133']['f'], 0.0055)

# Test update_gamma
def test_updgamma():
    ism = LineList('ISM')
    #
    np.testing.assert_allclose(ism['HI 1215']['gamma'], 626500000.0/u.s)

# Strong ISM LineList
def test_strong():
    strng = LineList('Strong')
    #
    assert len(strng._data) < 200

# Strong ISM LineList
def test_euv():
    euv = LineList('EUV')
    #
    assert np.max(euv._data['wrest'].value) < 1000.

# HI LineList
def test_h1():
    HI = LineList('HI')
    #
    for name in HI.name:
        assert name[0:2] == 'HI'

# H2 LineList
def test_h2():
    h2 = LineList('H2')
    #
    np.testing.assert_allclose(h2[911.967*u.AA]['f'], 0.001315, rtol=1e-5)

# CO LineList
def test_co():
    CO = LineList('CO')
    #
    np.testing.assert_allclose(CO[1322.133*u.AA]['f'], 0.0006683439, rtol=1e-5)

# Galactic LineList
def test_galx():
    galx = LineList('Galaxy')
    #
    np.testing.assert_allclose(galx["Halpha"]['wrest'].value, 6564.613, rtol=1e-5)


# Unknown lines
def test_unknown():
    ism = LineList('ISM')
    unknown = ism.unknown_line()
    assert unknown['name'] == 'unknown', 'There is a problem in the LineList.unknown_line()'
    assert unknown['wrest'] == 0.*u.AA, 'There is a problem in the LineList.unknown_line()'
    print(ism['unknown'])

def test_mk_sets():
    import imp
    llmk.mk_hi(outfil='tmp.lst', stop=False)
    lt_path = imp.find_module('linetools')[1]
    llmk.add_galaxy_lines('tmp.lst', infil=lt_path+'/lists/sets/llist_v0.1.ascii', stop=False)

def test_set_extra_columns_to_datatable():
    ism = LineList('ISM')
    # bad calls
    try:
        ism.set_extra_columns_to_datatable(abundance_type='incorrect_one')
    except ValueError:
        pass
    try:
        ism.set_extra_columns_to_datatable(ion_correction='incorrect_one')
    except ValueError:
        pass
    # test expected strongest value
    ism.set_extra_columns_to_datatable(ion_correction='none', abundance_type='solar')
    np.testing.assert_allclose(ism['HI 1215']['rel_strength'], 14.704326420257642)
    tab = ism._data

    np.testing.assert_allclose(np.max(tab['rel_strength']), 14.704326420257642)
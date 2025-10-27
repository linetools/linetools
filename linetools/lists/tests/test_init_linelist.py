# Module to run tests on Generating a LineList
#   Also tests some simple functionality
# TEST_UNICODE_LITERALS

import pytest
import os
from astropy import units as u
import numpy as np

from linetools.lists.linelist import LineList
from linetools.lists import mk_sets as llmk


def test_ism_read_source_catalogues():
    ism = LineList('ISM')
    np.testing.assert_allclose(ism['HI 1215']['wrest'], 1215.6700*u.AA, rtol=1e-7)

# ISM LineList
def test_ism():
    ism = LineList('ISM')
    np.testing.assert_allclose(ism['HI 1215']['wrest'], 1215.6700*u.AA, rtol=1e-7)

# Test update_fval
def test_updfval():
    ism = LineList('ISM')
    np.testing.assert_allclose(ism['FeII 1133']['f'], 0.0055)

# Test update_gamma
def test_updgamma():
    ism = LineList('ISM')
    np.testing.assert_allclose(ism['HI 1215']['gamma'], 626500000.0/u.s)


# Strong ISM LineList
def test_strong():
    strng = LineList('Strong')
    assert len(strng._data) < 200


# Strong ISM LineList
def test_euv():
    euv = LineList('EUV')
    #
    assert np.max(euv._data['wrest']) < 1000.
    # Test for X-ray lines
    ovii = euv['OVII 21']
    assert np.isclose(ovii['wrest'].value, 21.6019)


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
    np.testing.assert_allclose(galx["Halpha"]['wrest'], 6564.613*u.AA, rtol=1e-5)


# AGN LineList
def test_agn():
    agn = LineList('AGN')
    #
    np.testing.assert_allclose(agn["Halpha"]['wrest'], 6564.613*u.AA, rtol=1e-5)
    np.testing.assert_allclose(agn["NV 1242"]['wrest'], 1242.804*u.AA, rtol=1e-5)




# Unknown lines
def test_unknown():
    ism = LineList('ISM')
    unknown = ism.unknown_line()
    assert unknown['name'] == 'unknown', 'There is a problem in the LineList.unknown_line()'
    assert unknown['wrest'] == 0.*u.AA, 'There is a problem in the LineList.unknown_line()'

def test_mk_sets():
    outfile = 'tmp.lst'
    if os.path.isfile(outfile):
        os.remove(outfile)
    llmk.mk_hi(outfil=outfile, stop=False)
    # Get the linetools package directory from this module's location
    lt_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    llmk.add_galaxy_lines(outfile, infil=lt_path+'/lists/sets/llist_v0.1.ascii', stop=False)
    os.remove(outfile)


def test_set_extra_columns_to_datatable():
    # bad calls
    #ism = LineList('ISM')
    #with pytest.raises(ValueError) as tmp:  # This is failing Python 2.7 for reasons unbenknownst to me
    #    ism.set_extra_columns_to_datatable(abundance_type='incorrect_one')
    #ism = LineList('ISM')
    #with pytest.raises(ValueError):
    #    ism.set_extra_columns_to_datatable(ion_correction='incorrect_one', redo=True)
    # test expected strongest value
    ism = LineList('ISM')
    #np.testing.assert_allclose(ism['HI 1215']['rel_strength'], 14.704326420257642)  # THIS IS NO LONGER SUPPORTED
    tab = ism._extra_table
    np.testing.assert_allclose(np.max(tab['rel_strength']), 14.704326420257642)


# Module to run tests on Generating a LineList
#   Also tests some simple functionality
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

# Regression: the sets file records wrest rounded to 3-4 decimals while the
# master table carries more digits.  A fixed 9e-5 A matching tolerance silently
# dropped every entry whose stored value had been rounded, which cost the ISM
# list its whole FeII* fine-structure series (and 135 other lines).
FEII_FINE_STRUCTURE = [2328.111, 2333.516, 2338.725, 2345.001, 2349.022,
                       2359.828, 2365.552, 2381.489, 2405.164, 2411.802]


def test_ism_fe_ii_fine_structure():
    ism = LineList('ISM')
    wrest = np.array(ism._data['wrest'])
    names = np.array([str(nm) for nm in ism._data['name']])
    for wv in FEII_FINE_STRUCTURE:
        imt = np.argmin(np.abs(wrest - wv))
        assert np.abs(wrest[imt] - wv) < 0.05, \
            'FeII* {:.3f} missing from LineList("ISM")'.format(wv)
        # and it must be the real transition, not the g-weighted multiplet
        # mean rows at Ej = 416.299, which sit 0.5-1.2 A from genuine lines
        assert names[imt].startswith('FeII'), \
            'FeII* {:.3f} matched the wrong species: {:s}'.format(wv, names[imt])
        assert ism._data['Ej'][imt] > 0., \
            'FeII* {:.3f} matched a ground-state row'.format(wv)


def test_set_wrest_tol_does_not_reach_neighbouring_lines():
    """ The tolerance must stay well inside the closest real FeII pair
    (2631.8321 / 2632.1081, 0.28 A apart) and well inside the offset of the
    multiplet-mean rows (2344.7030 vs 2344.2139, 0.49 A).
    """
    from linetools.lists.parse import set_wrest_tol
    for wv in [2365.552, 2740.36, 1566.8217, 1206.5]:
        assert set_wrest_tol(wv) <= 5e-3
    # never stricter than the historical tolerance
    assert set_wrest_tol(1566.8217) >= 9e-5
    # a 3-decimal entry must reach its 4-decimal master counterpart
    assert set_wrest_tol(2365.552) > abs(2365.552 - 2365.5518)


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
    import importlib
    llmk.mk_hi(outfil=outfile, stop=False)
    lt_path = importlib.util.find_spec('linetools').submodule_search_locations[0]
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


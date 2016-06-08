# Module to run tests on RelAbund class

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pdb
import pytest
from astropy.table import Table
import numpy as np

from ..relabund import RelAbund

def make_class():
    XY = RelAbund()
    XY._data = {6: dict(flag=1, XH=-1., sigXH=0.2, sig=0.05),
                8: dict(flag=2, XH=-1.4, sigXH=0.25, sig=0.05),
                14: dict(flag=1, XH=-1.1, sigXH=0.25, sig=0.05),
                26: dict(flag=1, XH=-1.4, sigXH=0.25, sig=0.05),
                32: dict(flag=3, XH=-0.8, sigXH=0.25, sig=0.05),
                }
    return XY

def test_init():
    XY = RelAbund()
    assert XY.solar_ref == 'Asplund2009'
    assert isinstance(XY._data, dict)


def test_item():
    # Generate
    XY = make_class()
    # Test simple XH item
    CH = XY[6]
    np.testing.assert_allclose(CH['val'], -1.)
    # Test ratio
    SiFe = XY[14,26]
    np.testing.assert_allclose(SiFe['val'], +0.3)
    np.testing.assert_allclose(SiFe['sig'], 0.05*np.sqrt(2))
    # Test ratio with limits
    GeFe = XY['Ge','Fe']
    assert GeFe['flag'] == 3
    OFe = XY[8,26]
    assert OFe['flag'] == 2

def test_from_clmpair():
    CSi = RelAbund.from_clm_pair('C', 13.5, 'Si', 13.2)
    np.testing.assert_allclose(CSi[6,14]['val'], -0.62)

def test_from_iontbl():
    # Low-ions
    tbl = Table()
    tbl['Z'] = [6,6,8,14,26]
    tbl['ion'] = [2,4,1,2,2]
    tbl['flag_N'] = [2,1,2,1,1]
    tbl['logN'] = [15., 13., 15.5, 14., 13.]
    tbl['sig_logN'] = [0.05]*5
    # Instantiate without Ej
    XY = RelAbund.from_ionclm_table((1,20.5,0.2), tbl)
    # Test
    assert len(XY._data.keys()) == 4
    CH = XY[6]
    # Now with Ej
    tbl['Ej'] = [63., 0., 0., 0., 0.]
    XY = RelAbund.from_ionclm_table((1,20.5,0.2), tbl)
    assert len(XY._data.keys()) == 3
    with pytest.raises(KeyError):
        XY[6]

def test_table():
    # Generate
    XY = make_class()
    #
    tbl = XY.table()
    assert len(tbl) == 5

# Module to run tests on initializing AbsSurvey

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssystem import GenericAbsSystem
from linetools.isgm.abssurvey import GenericAbsSurvey

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_generic():
    # Start
    gensurvey = GenericAbsSurvey()

    # Systems
    coord = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys1 = GenericAbsSystem(coord, 1.244, [-300,300.]*u.km/u.s, NHI=16.)
    gensys1.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143*u.deg, dec=42.4321*u.deg)
    gensys2 = GenericAbsSystem(coord2, 1.744, [-300,300.]*u.km/u.s, NHI=17.)
    gensys2.name = 'Sys2'
    # Combine
    gensurvey.nsys = 2
    gensurvey._abs_sys.append(gensys1)
    gensurvey._abs_sys.append(gensys2)

    # Attribute
    aNHI = gensurvey.NHI
    np.testing.assert_allclose(aNHI, np.array([16.,17.]))

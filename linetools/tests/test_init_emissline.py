# Module to run tests on initializing EmLine

# TEST_UNICODE_LITERALS
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pytest
from astropy import units as u

from linetools.spectralline import EmLine, SpectralLine
from linetools import utils as ltu

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_mk_emissline():
    # Init Halpha
    emisslin = EmLine(6564.613*u.AA)
    assert emisslin.name == 'Halpha'

    # Init Halpha by name
    emisslin2 = EmLine('Halpha')
    np.testing.assert_allclose(emisslin2.data['wrest'], 6564.613*u.AA)


def test_dicts():
    # Init Halpha
    emisslin = EmLine(6564.613*u.AA)
    emisslin.analy['spec'] = 'tmp.fits'
    edict = emisslin.to_dict()
    assert isinstance(edict, dict)
    # Write
    ltu.savejson('tmp.json', edict, overwrite=True)
    # Read
    newdict = ltu.loadjson('tmp.json')
    newlin = SpectralLine.from_dict(newdict)
    assert newlin.name == 'Halpha'
    assert newlin.ltype == 'Em'

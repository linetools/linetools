# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from linetools.spectralline import AbsLine, SpectralLine
from linetools import utils as ltu

'''
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)
'''

def test_mk_absline():
    # Init HI Lya
    abslin = AbsLine(1215.6700*u.AA)
    np.testing.assert_allclose(abslin.data['f'], 0.4164)

    # Init CII 1334 with LineList
    abslin2 = AbsLine(1334.5323*u.AA, linelist='Strong')
    np.testing.assert_allclose(abslin2.data['Ek'], 74932.62 / u.cm)

    # Init CII 1334 by name
    abslin3 = AbsLine('CII 1334')
    np.testing.assert_allclose(abslin3.data['wrest'], 1334.5323*u.AA)


def test_dicts():
    # Init HI Lya
    abslin = AbsLine(1215.6700*u.AA)
    abslin.analy['spec'] = 'tmp.fits'
    adict = abslin.to_dict()
    assert isinstance(adict, dict)
    # Write
    #pdb.set_trace()
    ltu.savejson('tmp.json', adict, overwrite=True)
    # Read
    newdict = ltu.loadjson('tmp.json')
    newlin = SpectralLine.from_dict(newdict)
    assert newlin.name == 'HI 1215'

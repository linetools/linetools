# Module to run tests on Generating a LineList
#   Also tests some simple functionality

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import os
import pdb
import pytest
from astropy import units as u

from linetools.guis import utils as ltgu
from linetools.spectra import io as lsio


# Set of Input lines
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../../spectra/tests/files')
    return os.path.join(data_dir, filename)


# ISM LineList
def test_navigate():
    # Init
    nav_dict = dict(nav=ltgu.navigate(0,0,init=True))
    assert isinstance(nav_dict['nav'], list)
    nav_dict['xmnx'] = [0., 1]
    nav_dict['ymnx'] = [0., 1]
    nav_dict['sv_xy'] = [[0,1], [0,1]]
    nav_dict['tmp_xy'] = None
    # Usage
    o = type(str('Dummy'), (object,), {})
    o.xdata = 22.
    o.ydata = 1.
    for key in nav_dict['nav']:
        o.key = key
        ltgu.navigate(nav_dict, o)


def test_doublet():
    o = type(str('Dummy'), (object,), {})
    o.xdata = 5000.
    i = type(str('Dummy2'), (object,), {})
    for key in ['C','M','4','X','8','B']:
        o.key = key
        _ = ltgu.set_doublet(i, o)


def test_llist():
    # Init
    idict = ltgu.set_llist('Strong')
    idict = ltgu.set_llist([1215.670*u.AA])
    assert idict['List'] == 'input.lst'


def test_rdspec():
    spec, spec_fil = ltgu.read_spec(data_path('UM184_nF.fits'))
    #
    ispec = lsio.readspec(data_path('UM184_nF.fits'))
    spec, spec_fil = ltgu.read_spec(ispec)


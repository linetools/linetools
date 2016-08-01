# Module to run tests on Generating a LineList
#   Also tests some simple functionality

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import os, sys
import pdb
import pytest
from astropy import units as u

from PyQt4 import QtGui

from linetools.guis import xspecgui, xabssysgui
from linetools.guis import utils as ltgu
from linetools.spectra import io as lsio
from linetools.isgm.abssystem import GenericAbsSystem

app = QtGui.QApplication(sys.argv)

# Set of Input lines
def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../../spectra/tests/files')
    return os.path.join(data_dir, filename)


def test_navigate():
    # Init
    nav_dict = dict(nav=ltgu.navigate(0,0,init=True))
    assert isinstance(nav_dict['nav'], list)
    nav_dict['x_minmax'] = [0., 1]
    nav_dict['y_minmax'] = [0., 1]
    nav_dict['sv_xy_minmax'] = [[0,1], [0,1]]
    nav_dict['tmp_xy'] = None
    # Usage
    o = type(str('Dummy'), (object,), {})
    o.xdata = 22.
    o.ydata = 1.
    for key in nav_dict['nav']:
        o.key = key
        ltgu.navigate(nav_dict, o)
    # test wrong key event
    o.xdata = 'this_is_not_float'
    out = ltgu.navigate(nav_dict, o)
    assert out == 0
    # test event 's'
    o.key = 's'
    nav_dict['tmp_xy'] = (22, 1) #  i.e. not None
    ltgu.navigate(nav_dict, o)
    # test event 'y'
    o.key = 'y'
    ltgu.navigate(nav_dict, o, wave = np.linspace(1000,2000,100), flux = np.ones(100))


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
    idict = ltgu.set_llist('None')
    idict = ltgu.set_llist('OVI')
    # wrong format
    with pytest.raises(IOError):
        idict = ltgu.set_llist((1,2))  # input is a tuple, so it is wrong.


def test_rdspec():
    spec, spec_fil = ltgu.read_spec(data_path('UM184_nF.fits'))
    #
    ispec = lsio.readspec(data_path('UM184_nF.fits'))
    spec, spec_fil = ltgu.read_spec(ispec)


def test_xspecgui():
    # Init
    spec_fil = data_path('UM184_nF.fits')
    xsgui = xspecgui.XSpecGui(spec_fil, unit_test=True)


def test_xabsgui():
    # Init
    spec_fil = data_path('UM184_nF.fits')
    abs_sys = GenericAbsSystem((0.,0.), 3., [-500,500]*u.km/u.s)
    xabsgui = xabssysgui.XAbsSysGui(spec_fil, abs_sys)

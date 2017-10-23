# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.isgm.io import read_joebvp_to_components
from linetools.isgm.io import write_joebvp_from_components

from .utils import mk_comp, compare_two_files

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_read_joebvp():
    vp_file = data_path('group_9.VP')
    # Load
    icoord = SkyCoord(ra=12., dec=-12, unit='deg')
    comps = read_joebvp_to_components(vp_file, icoord)
    # Test
    assert isinstance(comps, list)
    assert len(comps) == 2
    assert len(comps[1]._abslines) == 2
    assert comps[1]._abslines[0].attrib['N'].value > 0.


def test_complist_to_joebvp():
    # will write a file in directory ./files/
    abscomp, HIlines = mk_comp('HI', b=15*u.km/u.s, use_rand=False)
    comp_list = [abscomp, abscomp]
    write_joebvp_from_components(comp_list, 'test.fits', data_path('test_joebvp_repr.joebvp'))
    # now read the output and compare to reference
    compare_two_files(data_path('test_joebvp_repr.joebvp'),
                      resource_filename('linetools', '/data/tests/test_joebvp_repr_reference.joebvp'))
    # now add attribute to comp and compare again
    abscomp.attrib['b'] = 15*u.km/u.s
    write_joebvp_from_components(comp_list, 'test.fits', data_path('test_joebvp_repr.joebvp'))
    compare_two_files(data_path('test_joebvp_repr.joebvp'),
                      resource_filename('linetools', '/data/tests/test_joebvp_repr_reference.joebvp'))


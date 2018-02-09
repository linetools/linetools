# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os

from pkg_resources import resource_filename

from astropy.coordinates import SkyCoord
from astropy import units as u

import linetools.isgm.io as ltiio

from linetools.isgm.tests.utils import mk_comp, compare_two_files

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_read_joebvp():
    vp_file = data_path('group_9.VP')
    # Load
    icoord = SkyCoord(ra=12., dec=-12, unit='deg')
    comps = ltiio.read_joebvp_to_components(vp_file, icoord)
    # Test
    assert isinstance(comps, list)
    assert len(comps) == 2
    assert len(comps[1]._abslines) == 2
    assert comps[1]._abslines[0].attrib['N'].value > 0.


def test_complist_to_joebvp():
    # will write a file in directory ./files/
    abscomp, HIlines = mk_comp('HI', b=15*u.km/u.s, use_rand=False)
    comp_list = [abscomp, abscomp]
    ltiio.write_joebvp_from_components(comp_list, 'test.fits', data_path('test_joebvp_repr.joebvp'))
    # now read the output and compare to reference
    compare_two_files(data_path('test_joebvp_repr.joebvp'),
                      resource_filename('linetools', '/data/tests/test_joebvp_repr_reference.joebvp'))
    # now add attribute to comp and compare again
    abscomp.attrib['b'] = 15*u.km/u.s
    ltiio.write_joebvp_from_components(comp_list, 'test.fits', data_path('test_joebvp_repr.joebvp'))
    compare_two_files(data_path('test_joebvp_repr.joebvp'),
                      resource_filename('linetools', '/data/tests/test_joebvp_repr_reference.joebvp'))


# def test_read_write_igmg():
if 1:
    # read
    igmg_file = data_path('J1410+2304_model.json')
    comps = ltiio.read_igmg_to_components(igmg_file)
    assert comps[0].name == 'OIII_z0.53512'
    assert comps[0].reliability == 'a'
    assert comps[8].zcomp == 0.068
    assert len(comps) == 132
    # write
    ltiio.write_igmg_from_components(comps[:2], specfile='test.fits', fwhm=3, outfile=data_path('IGM_model.json'))
    compare_two_files(data_path('IGM_model.json'),
                      resource_filename('linetools', '/data/tests/IGM_model_reference.json'), except_l2_has='2018-Feb-09')

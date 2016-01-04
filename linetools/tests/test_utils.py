from __future__ import print_function, absolute_import, division, unicode_literals

from ..utils import between, v_from_z, savejson, loadjson, z_from_v
from ..utils import radec_to_coord
import numpy as np
from astropy import units as u
import pdb

def test_between():
    x = [1,2,3,4,5]
    c0 = between(x, 2.5, 3.5)
    assert np.all(c0 == np.array([False, False, True, False, False]))

def test_zfromv():
    z = z_from_v(2., 1000.)
    #
    np.testing.assert_allclose(z, 2.0100236684175417)

def test_vfromz():
    v = v_from_z(2., 2.1)
    #
    assert v.unit == u.km/u.s
    np.testing.assert_allclose(v.value, -9826.62006340679)

def test_save_load_json():
    tmp_dict = dict(a=1, b=2, c='adsf')
    # Write
    savejson('tmp.json', tmp_dict, overwrite=True)
    # Load
    new_dict = loadjson('tmp.json')
    assert new_dict['a'] == 1


def test_radeccoord():
    for radec in ['J124511+144523', '124511+144523',
                  'J12:45:11+14:45:23', ('12:45:11','+14:45:23'),
                  ('12:45:11', '14:45:23'), ('12 45 11','+14 45 23')]:
        coord = radec_to_coord(radec)
        # Test
        np.testing.assert_allclose(coord.ra.value, 191.2958333333333)

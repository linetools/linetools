from __future__ import print_function, absolute_import, division, unicode_literals

from ..utils import between, v_from_z, savejson, loadjson, z_from_v
from ..utils import radec_to_coord, convert_quantity_in_dict
import numpy as np
from astropy import units as u
import linetools.utils as ltu
import pdb

def test_convert_qdict():
    vlim = dict(unit='km/s', value=50.)
    time = dict(unit='s', value=100.)
    mass = dict(unit='g', value=10.)
    idict = dict(vlim=vlim, nest=dict(time=time, mass=mass))
    newdict = convert_quantity_in_dict(idict)
    assert newdict['nest']['mass'] == 10*u.g
    # Simple one
    obj = convert_quantity_in_dict(vlim)
    assert obj.unit == u.km/u.s

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
    # Write with gzip
    savejson('tmp.json.gz', tmp_dict, overwrite=True)
    # Load
    new_dict = loadjson('tmp.json.gz')
    assert new_dict['a'] == 1


def test_radeccoord():
    for radec in ['J124511+144523', '124511+144523',
                  'J12:45:11+14:45:23', ('12:45:11','+14:45:23'),
                  ('12:45:11', '14:45:23'), ('12 45 11','+14 45 23')]:
        coord = radec_to_coord(radec)
        # Test
        np.testing.assert_allclose(coord.ra.value, 191.2958333333333)


def test_overlapping_chunks():
    chunk1 = (1,2,3,4)
    chunk2 = [3,4,5,6]
    t = ltu.overlapping_chunks(chunk1, chunk2)
    assert t
    t = ltu.overlapping_chunks(chunk2, chunk1)
    assert t

    chunk2 = np.array([5,7])
    f = ltu.overlapping_chunks(chunk2, chunk1)
    assert ~f
    f = ltu.overlapping_chunks(chunk2*u.AA, chunk1*u.AA)
    assert ~f

    # Wrong format
    try:
        f = ltu.overlapping_chunks(chunk2*u.AA, chunk1)
    except ValueError:
        pass
    try:
        f = ltu.overlapping_chunks(chunk2*u.AA, chunk1*u.K)
    except ValueError:
        pass

    # not sorted
    try:
        f = ltu.overlapping_chunks((1,2,3,4), (4,5,3,6))
    except ValueError:
        pass

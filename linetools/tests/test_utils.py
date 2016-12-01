# Tests of linetools.utils

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import linetools.utils as ltu


def test_name_from_coord():
    coord = SkyCoord(ra=15.4, dec=-23.1, unit='deg')
    name = ltu.name_from_coord(coord)
    assert name == 'J010136.00-230600.0'
    # Name 2
    name2 = ltu.name_from_coord(coord, precision=(0, 0))
    assert name2 == 'J010136-230600'


def test_convert_qdict():
    vlim = dict(unit='km/s', value=50.)
    time = dict(unit='s', value=100.)
    mass = dict(unit='g', value=10.)
    idict = dict(vlim=vlim, nest=dict(time=time, mass=mass))
    newdict = ltu.convert_quantity_in_dict(idict)
    assert newdict['nest']['mass'] == 10 * u.g
    # Simple one
    obj = ltu.convert_quantity_in_dict(vlim)
    assert obj.unit == u.km / u.s


def test_between():
    x = [1, 2, 3, 4, 5]
    c0 = ltu.between(x, 2.5, 3.5)
    assert np.all(c0 == np.array([False, False, True, False, False]))


def test_dz_from_dv():
    dz = ltu.dz_from_dv(1000. * u.km / u.s, 2.)
    np.testing.assert_allclose(dz, 0.0100236684175417)
    dz1 = ltu.dz_from_dv([1000., 1000., 1000.] * u.km / u.s, np.array([2., 2., 2.]))
    dz2 = ltu.dz_from_dv([1000., 1000., 1000.] * u.km / u.s, [2., 2., 2.])
    np.testing.assert_allclose(dz1, [0.0100236684175417] * 3)
    np.testing.assert_allclose(dz2, [0.0100236684175417] * 3)

    # non-relativistic
    dz = ltu.dz_from_dv(1000. * u.km / u.s, 2., rel=False)
    np.testing.assert_allclose(dz, 0.010006922855944561)

    # test expected errors
    with pytest.raises(IOError):
        ltu.dz_from_dv('dv_not_a_quantity', 2.)
    with pytest.raises(IOError):
        ltu.dz_from_dv(1000. * u.km / u.s, 'zref_not_a_float_nor_array')
    with pytest.raises(IOError):
        ltu.dz_from_dv([1000., 1000., 1000.] * u.km / u.s, np.array([2., 2.]))  # wrong shape for zref
    with pytest.raises(IOError):
        ltu.dz_from_dv([1000., 1000., 1000.] * u.km, 1.)  # wrong dv units


def test_z_from_dv():
    z = ltu.z_from_dv(1000. * u.km / u.s, 2.)
    np.testing.assert_allclose(z, 2.0100236684175417)


def test_dv_from_z():
    dv = ltu.dv_from_z(2.1, 2.)
    assert dv.unit == u.km / u.s
    np.testing.assert_allclose(dv.value, 9826.620063406788, rtol=1e-6)
    dv = ltu.dv_from_z([2.1, 2.1, 2.1], 2.)
    np.testing.assert_allclose(dv.value, [9826.620063406788] * 3, rtol=1e-6)
    # non-relativistic
    dv = ltu.dv_from_z(2.1, 2., rel=False)
    np.testing.assert_allclose(dv.value, 9993.08193333334, rtol=1e-6)

    # test expected errors
    with pytest.raises(IOError):
        ltu.dv_from_z('z_not_a_float_or_array', 1.)
    with pytest.raises(IOError):
        ltu.dv_from_z(2.1, 'zref_not_a_float_nor_array')
    with pytest.raises(IOError):
        ltu.dv_from_z(np.array([2.1, 2.1, 2.1]), np.array([2., 2.]))  # wrong shape for zref


def test_save_load_json():
    tmp_dict = dict(a=1, b=2, c='adsf')
    # Write
    ltu.savejson('tmp.json', tmp_dict, overwrite=True)
    # Load
    new_dict = ltu.loadjson('tmp.json')
    assert new_dict['a'] == 1
    # Write with gzip
    ltu.savejson('tmp.json.gz', tmp_dict, overwrite=True)
    # Load
    new_dict = ltu.loadjson('tmp.json.gz')
    assert new_dict['a'] == 1
    # Write with easy to read
    ltu.savejson('tmp2.json', tmp_dict, overwrite=True, easy_to_read=True)
    new_dict2 = ltu.loadjson('tmp2.json')
    assert new_dict2['a'] == 1


def test_radeccoord():
    for radec in ['J124511+144523', '124511+144523',
                  'J12:45:11+14:45:23', ('12:45:11', '+14:45:23'),
                  ('12:45:11', '14:45:23'), ('12 45 11', '+14 45 23')]:
        coord = ltu.radec_to_coord(radec)
        # Test
        np.testing.assert_allclose(coord.ra.value, 191.2958333333333)


def test_overlapping_chunks():
    chunk1 = (1, 2, 3, 4)
    chunk2 = [3, 4, 5, 6]
    t = ltu.overlapping_chunks(chunk1, chunk2)
    assert t
    t = ltu.overlapping_chunks(chunk2, chunk1)
    assert t

    chunk2 = np.array([5, 7])
    f = ltu.overlapping_chunks(chunk2, chunk1)
    assert ~f
    f = ltu.overlapping_chunks(chunk2 * u.AA, chunk1 * u.AA)
    assert ~f

    # Wrong format
    try:
        f = ltu.overlapping_chunks(chunk2 * u.AA, chunk1)
    except ValueError:
        pass
    try:
        f = ltu.overlapping_chunks(chunk2 * u.AA, chunk1 * u.K)
    except ValueError:
        pass

    # not sorted
    try:
        f = ltu.overlapping_chunks((1, 2, 3, 4), (4, 5, 3, 6))
    except ValueError:
        pass

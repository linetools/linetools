# Module to run tests on generating AbsSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import AbsSystem, GenericAbsSystem, LymanAbsSystem
from linetools.spectralline import AbsLine
from linetools.spectra import io

import pdb

try:
    unicode
except NameError:
    unicode = str


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../../spectra/tests/files')
    return os.path.join(data_dir, filename)


def init_system():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lya.attrib['coord'] = radec
    lyb = AbsLine(1025.7222*u.AA, z=lya.z)
    lyb.limits.set([-300.,300.]*u.km/u.s)
    lyb.attrib['coord'] = radec
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans, z=2.92939)
        iline.attrib['coord'] = radec
        iline.limits.set([-250.,80.]*u.km/u.s)
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines)
    # Instantiate
    gensys = GenericAbsSystem.from_components([abscomp,SiII_comp])
    return gensys


def test_list_of_abslines():
    spec = io.readspec(data_path('UM184_nF.fits'))

    gensys = init_system()
    # Now the list
    abslines = gensys.list_of_abslines()
    # Test
    assert len(abslines) == 6
    # Measure AODM
    gensys.measure_aodm(spec=spec)
    # trans table
    gensys.fill_trans()
    assert len(gensys._trans) == 6
    # Grab one line
    lyb = gensys.get_absline('HI 1025')
    np.testing.assert_allclose(lyb.wrest.value, 1025.7222)
    lyb = gensys.get_absline(1025.72*u.AA)
    np.testing.assert_allclose(lyb.wrest.value, 1025.7222)

def test_measure_ew():
    spec = io.readspec(data_path('UM184_nF.fits'))
    gensys = init_system()
    # Measure EWs
    gensys.measure_restew(spec=spec)
    for aline in gensys.list_of_abslines():
        assert aline.attrib['flag_EW'] == 1
    # Skip do_analysis = 0
    gensys2 = init_system()
    lyb = gensys2.get_absline('HI 1025')
    lyb.analy['do_analysis'] = 0
    gensys2.measure_restew(spec=spec)
    assert lyb.attrib['flag_EW'] == 0

def test_ionn():
    gensys = init_system()
    # Component columns
    gensys.update_component_colm()
    # ionN
    gensys.fill_ionN()
    assert len(gensys._ionN) == 2
    #
    gensys.NHI = 15.3
    gensys.sig_NHI = 0.3
    gensys.flag_NHI = 1
    gensys.fill_ionN(NHI_obj=gensys)
    HI = np.where(gensys._ionN['Z']==1)
    np.testing.assert_allclose(gensys._ionN[HI]['logN'], 15.3)


def test_get_component():
    gensys = init_system()
    # Grab SiII
    SiII = gensys.get_component((14,2))
    assert isinstance(SiII, AbsComponent)
    # Fail
    SiV = gensys.get_component((14,5))
    assert SiV is None
    # Line
    alines = gensys.list_of_abslines()
    comp = gensys.get_component(alines[0])
    assert isinstance(comp, AbsComponent)
    # Junk
    with pytest.raises(IOError):
        junk = gensys.get_component(None)


def test_todict():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lya.attrib['coord'] = radec
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    lyb.attrib['coord'] = radec
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # Instantiate
    HIsys = LymanAbsSystem.from_components([abscomp])
    # Dict
    adict = HIsys.to_dict()
    assert isinstance(adict, dict)
    # Instantiate
    #pdb.set_trace()
    newsys = AbsSystem.from_dict(adict)
    assert isinstance(newsys, AbsSystem)


@pytest.mark.skipif("sys.version_info >= (3,0)")
def test_todict_withjson():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lya.attrib['z'] = 2.92939
    lya.attrib['coord'] = radec
    lyb = AbsLine(1025.7222*u.AA)
    lyb.limits.set([-300.,300.]*u.km/u.s)
    lyb.attrib['z'] = lya.attrib['z']
    lyb.attrib['coord'] = radec
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # Instantiate
    HIsys = LymanAbsSystem.from_components([abscomp])
    # Dict
    adict = HIsys.to_dict()
    assert isinstance(adict, dict)
    # Verify it is JSON compatible (failing in Python 3)
    import io,json
    with io.open('tmp.json', 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(adict, sort_keys=True, indent=4,
                                   separators=(',', ': '))))

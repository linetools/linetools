# Module to run tests on generating AbsSystem

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

from linetools.isgm.abscomponent import AbsComponent
from linetools.isgm.abssystem import AbsSystem, GenericAbsSystem, LymanAbsSystem
from linetools.spectralline import AbsLine

import pdb

def test_list_of_abslines():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines)
    SiII_comp.coord = radec
    # Instantiate
    gensys = GenericAbsSystem.from_components([abscomp,SiII_comp])
    # Now the list
    abslines = gensys.list_of_abslines()
    # Test
    assert len(abslines) == 6
    # Grab one line
    lyb = gensys.absline('HI 1025')
    np.testing.assert_allclose(lyb.wrest.value, 1025.7222)
    lyb = gensys.absline(1025.72*u.AA)
    np.testing.assert_allclose(lyb.wrest.value, 1025.7222)


def test_todict():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # Instantiate
    HIsys = LymanAbsSystem.from_components([abscomp])
    # Dict
    adict = HIsys.to_dict()
    assert isinstance(adict, dict)
    # Verify it is JSON compatible
    import io,json
    with io.open('tmp.json', 'w', encoding='utf-8') as f:
        f.write(unicode(json.dumps(adict, sort_keys=True, indent=4,
                                   separators=(',', ': '))))
    # Instantiate
    newsys = AbsSystem.from_dict(adict)
    assert isinstance(newsys, AbsSystem)

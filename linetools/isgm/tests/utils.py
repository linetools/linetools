# Module to helping to run tests

from __future__ import print_function, absolute_import, division, unicode_literals

import os
import pdb

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssightline import GenericAbsSightline
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def lyman_comp(radec, z=2.92939):
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA, z=z)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lya.attrib['flag_N'] = 1
    lya.attrib['N'] = 1e17 /  u.cm**2
    lya.attrib['coord'] = radec
    lyb = AbsLine(1025.7222*u.AA, z=z)
    lyb.limits.set([-300.,300.]*u.km/u.s)
    lyb.attrib['coord'] = radec
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.synthesize_colm()

    return abscomp


def si2_comp(radec):
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
    SiII_comp.logN = 15.
    SiII_comp.flag_N = 1
    #
    return SiII_comp


def make_gensl():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    gensl = GenericAbsSightline.from_components([abscomp, SiII_comp])
    return gensl


def write_comps_to_sys():
    from linetools.isgm.abssystem import GenericAbsSystem
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec, z=2.92940)
    # SiII
    SiII_comp = si2_comp(radec)
    gensl = GenericAbsSystem.from_components([abscomp, SiII_comp])
    # Write
    gensl.write_json()

# Command line execution
if __name__ == '__main__':

    flg_tab = 0
    flg_tab += 2**0  # Comps to System

    # Generate
    if flg_tab & (2**0):
        write_comps_to_sys()
        # WRites to Foo.json



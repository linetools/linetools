# Module to helping to run tests

from __future__ import print_function, absolute_import, division, unicode_literals

import os
import numpy as np
import pdb

from pkg_resources import resource_filename

from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools.isgm.abssystem import GenericAbsSystem
from linetools.isgm.abssightline import GenericAbsSightline
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.spectra import io as lsio
from linetools.lists.linelist import LineList

ism = LineList('ISM')


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)




def lyman_comp(radec, z=2.92939):
    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA, z=z, linelist=ism)
    lya.limits.set([-300.,300.]*u.km/u.s)
    lya.attrib['flag_N'] = 1
    lya.attrib['N'] = 1e17 /  u.cm**2
    lya.attrib['sig_N'] = 1e16 /  u.cm**2
    lya.attrib['coord'] = radec
    # Lyb
    lyb = AbsLine(1025.7222*u.AA, z=z, linelist=ism)
    lyb.limits.set([-300.,300.]*u.km/u.s)
    lyb.attrib['coord'] = radec
    lyb.attrib['flag_N'] = 1
    lyb.attrib['N'] = 1e17 /  u.cm**2
    lyb.attrib['sig_N'] = 1e16 /  u.cm**2
    # Build
    abscomp = AbsComponent.from_abslines([lya,lyb])
    #abscomp.synthesize_colm()

    return abscomp


def si2_comp(radec, z=2.92939):
    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans, z=z, linelist=ism)
        iline.attrib['coord'] = radec
        iline.limits.set([-250.,80.]*u.km/u.s)
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines, skip_synth=True)
    SiII_comp.logN = 15.
    SiII_comp.flag_N = 1
    #
    return SiII_comp

def oi_comp(radec, vlim=[-250.,80.]*u.km/u.s, z=2.92939):
    # SiII
    OItrans = ['OI 1302']
    abslines = []
    for trans in OItrans:
        iline = AbsLine(trans, z=z, linelist=ism)
        iline.attrib['coord'] = radec
        iline.limits.set(vlim)
        abslines.append(iline)
    #
    OI_comp = AbsComponent.from_abslines(abslines, skip_synth=True)
    OI_comp.logN = 15.
    OI_comp.flag_N = 1
    #
    return OI_comp

def make_gensl():
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    # HI
    abscomp = lyman_comp(radec)
    # SiII
    SiII_comp = si2_comp(radec)
    gensl = GenericAbsSightline.from_components([abscomp, SiII_comp])
    # Add a system
    sys = GenericAbsSystem.from_components([abscomp, SiII_comp])
    gensl._abssystems = [sys]
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


def mk_comp(ctype,vlim=[-300.,300]*u.km/u.s,add_spec=False, use_rand=True,
            add_trans=False, zcomp=2.92939, b=20*u.km/u.s, **kwargs):
    # Read a spectrum Spec
    if add_spec:
        spec_file = resource_filename('linetools','/spectra/tests/files/UM184_nF.fits')
        xspec = lsio.readspec(spec_file)
    else:
        xspec = None
    # AbsLines
    if ctype == 'HI':
        all_trans = ['HI 1215', 'HI 1025']
    elif ctype == 'SiII':
        all_trans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
        if add_trans:
            all_trans += ['SiII 1193']
    elif ctype == 'SiII*':
        all_trans = ['SiII* 1264', 'SiII* 1533']
    elif ctype == 'SiIII':
        all_trans = ['SiIII 1206']
    abslines = []
    for trans in all_trans:
        iline = AbsLine(trans, z=zcomp, linelist=ism)
        if use_rand:
            rnd = np.random.rand()
        else:
            rnd = 0.
        iline.attrib['logN'] = 13.3 + rnd
        iline.attrib['sig_logN'] = 0.15
        iline.attrib['flag_N'] = 1
        iline.attrib['b'] = b
        iline.analy['spec'] = xspec
        iline.limits.set(vlim)
        _,_ = ltaa.linear_clm(iline.attrib)  # Loads N, sig_N
        abslines.append(iline)
    # Component
    abscomp = AbsComponent.from_abslines(abslines, **kwargs)
    return abscomp, abslines


def mk_comptable():
    from astropy.table import Table
    tab = Table()
    tab['ion_name'] = ['HI', 'HI', 'CIV', 'SiII', 'OVI']
    tab['Z'] = [1,1,4,14,8]
    tab['ion'] = [1,1,4,2,6]
    tab['z_comp'] = [0.05, 0.0999, 0.1, 0.1001, 0.6]
    tab['RA'] = [100.0] * len(tab) * u.deg
    tab['DEC'] = [-0.8] * len(tab) * u.deg
    tab['vmin'] = [-50.] * len(tab) * u.km/u.s
    tab['vmax'] = [100.] * len(tab) * u.km/u.s
    tab['Ej'] = [0.] *len(tab) / u.cm
    return tab


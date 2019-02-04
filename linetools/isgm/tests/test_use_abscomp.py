# Module to run tests on using AbsComponent

from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
import pdb

#import matplotlib
#matplotlib.use('Agg')

from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.isgm import utils as ltiu

from linetools.isgm.tests.utils import mk_comp, mk_comptable, ism

import imp, os
lt_path = imp.find_module('linetools')[1]

#import pdb
#pdb.set_trace()
# Set of Input lines

class DummyObj(object):
    pass


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_rest_limits():
    SiIItrans = ['SiII 1260', 'SiII 1304']
    radec = SkyCoord(ra=1., dec=1., unit='deg')
    abslines = []
    for kk,trans in enumerate(SiIItrans):
        iline = AbsLine(trans, z=1., linelist=ism)
        iline.attrib['coord'] = radec
        if kk == 0:
            iline.limits.set([-250.,80.]*u.km/u.s)
        else:
            iline.limits.set([-50.,180.]*u.km/u.s)
        abslines.append(iline)
    #
    comp = AbsComponent.from_abslines(abslines, chk_vel=False)
    comp.reset_limits_from_abslines()
    assert np.isclose(comp.vlim[1].value, 180.)

def test_write():
    # Component
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s)
    # Write
    abscomp.write(data_path('tmp.json'))
    # Read
    tmpcomp = abscomp.from_json(data_path('tmp.json'))
    assert tmpcomp.flag_N == abscomp.flag_N
    assert np.isclose(tmpcomp.limits.vmin.value, abscomp.limits.vmin.value)


def test_get_components_at_z():
    tab = mk_comptable()
    complist = ltiu.complist_from_table(tab)
    z01_comps = ltiu.get_components_at_z(complist, 0.1, [-1000,1000]*u.km/u.s)
    assert len(z01_comps) == 3
    # check expected errors
    with pytest.raises(IOError):
        ltiu.get_components_at_z([1,2,3], 0.1, [-1000,1000]*u.km/u.s) # wrong complist
    with pytest.raises(IOError):
        ltiu.get_components_at_z(complist, 0.1, [-1000,1000, 1000]*u.km/u.s) # wrong vlims size
    with pytest.raises(IOError):
        ltiu.get_components_at_z(complist, 0.1, [-1000,1000]*u.km)  # wrong vlims units


def test_cog():
    # Component
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # COG
    COG_dict = abscomp.cog(redo_EW=True)
    # Test
    np.testing.assert_allclose(COG_dict['logN'],13.693355878125537)
    np.testing.assert_allclose(COG_dict['sig_logN'],0.054323725737309987)



"""
def test_stack_plot():
    # AbsLine(s)
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    # Spectra
    xspec = lsio.readspec(lt_path+'/spectra/tests/files/UM184_nF.fits')
    lya.analy['spec'] = xspec
    lyb.analy['spec'] = xspec
    # Instantiate
    abscomp = AbsComponent.from_abslines([lya,lyb])
    # Plot
    abscomp.stack_plot(show=False)
"""


def test_repr_vpfit():
    abscomp, HIlines = mk_comp('HI', use_rand=False)
    s = abscomp.repr_vpfit()
    assert s == 'HI 2.92939 0.00000 10.00 0.00 13.30 0.00\n'

    s = abscomp.repr_vpfit(b=15*u.km/u.s)
    assert s == 'HI 2.92939 0.00000 15.00 0.00 13.30 0.00\n'

    abscomp.comment = 'Something'
    s = abscomp.repr_vpfit()
    assert s == 'HI 2.92939 0.00000 10.00 0.00 13.30 0.00! Something\n'
    s = abscomp.repr_vpfit(tie_strs=('a', 'b', 'CD'), fix_strs=('', 'f', ''))
    assert s == 'HI 2.92939a 0.00000 10.00F 0.00 13.30cd 0.00! Something\n'

    abscomp, SiIIlines = mk_comp('SiII', use_rand=False)
    s = abscomp.repr_vpfit()
    assert s == 'SiII 2.92939 0.00000 10.00 0.00 13.30 0.00\n'

    # errors
    with pytest.raises(TypeError):
        s = abscomp.repr_vpfit(tie_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_vpfit(fix_strs='bad_format')
    with pytest.raises(SyntaxError):
        s = abscomp.repr_vpfit(fix_strs=('1','2','3','4','5'))


def test_repr_alis():
    abscomp, HIlines = mk_comp('HI', use_rand=False)
    s = abscomp.repr_alis()
    assert s == 'voigt   ion=1H_I 13.30 redshift=2.92939 0.0 1.0E+04\n'

    abscomp, SiIIlines = mk_comp('SiII', use_rand=False)
    s = abscomp.repr_alis(T_kin=10**5*u.K)
    assert s == 'voigt   ion=28Si_II 13.30 redshift=2.92939 0.0 1.0E+05\n'

    abscomp.comment = 'Something'
    s = abscomp.repr_alis()
    assert s == 'voigt   ion=28Si_II 13.30 redshift=2.92939 0.0 1.0E+04# Something\n'
    s = abscomp.repr_alis(tie_strs=('a', 'b', 'CD',''), fix_strs=('', 'f', '', ''))
    assert s == 'voigt   ion=28Si_II 13.30a redshift=2.92939F 0.0cd 1.0E+04# Something\n'

    # errors
    with pytest.raises(TypeError):
        s = abscomp.repr_alis(tie_strs='bad_format')
    with pytest.raises(TypeError):
        s = abscomp.repr_alis(fix_strs='bad_format')
    with pytest.raises(SyntaxError):
        s = abscomp.repr_alis(fix_strs=('1','2','3','4','5'))


def test_repr_joebvp():
    # test with b=0, should be replaced by b_default
    abscomp, HIlines = mk_comp('HI', b=0*u.km/u.s, use_rand=False)
    s = abscomp.repr_joebvp('test.fits', b_default=3.3*u.km/u.s)
    assert s == 'test.fits|1215.67000|2.92939000|13.3000|3.3000|0.0000|2|2|2|-300.0000|300.0000|4772.06378|4781.62408|2.92939000|HI|none|\n' \
                'test.fits|1025.72220|2.92939000|13.3000|3.3000|0.0000|2|2|2|-300.0000|300.0000|4026.43132|4034.49783|2.92939000|HI|none|\n'
    # test with b != 0
    abscomp, HIlines = mk_comp('HI', b=15*u.km/u.s, use_rand=False)
    s = abscomp.repr_joebvp('test.fits', b_default=3.3*u.km/u.s)
    assert s == 'test.fits|1215.67000|2.92939000|13.3000|15.0000|0.0000|2|2|2|-300.0000|300.0000|4772.06378|4781.62408|2.92939000|HI|none|\n' \
                'test.fits|1025.72220|2.92939000|13.3000|15.0000|0.0000|2|2|2|-300.0000|300.0000|4026.43132|4034.49783|2.92939000|HI|none|\n'
    # test with comment
    abscomp.comment = 'Something'
    s = abscomp.repr_joebvp('test.fits')
    assert s == 'test.fits|1215.67000|2.92939000|13.3000|15.0000|0.0000|2|2|2|-300.0000|300.0000|4772.06378|4781.62408|2.92939000|HI|none|Something\n' \
                'test.fits|1025.72220|2.92939000|13.3000|15.0000|0.0000|2|2|2|-300.0000|300.0000|4026.43132|4034.49783|2.92939000|HI|none|Something\n'



def test_get_wvobs_chunks():
    abscomp, HIlines = mk_comp('HI', zcomp=0., vlim=[0,10]*u.km/u.s)
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    np.testing.assert_allclose(wvobs_chunks[0][0], 1215.67*u.AA)
    np.testing.assert_allclose(wvobs_chunks[0][1], 1215.71055106*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][0], 1025.7222*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][1], 1025.75641498*u.AA)
    abscomp, HIlines = mk_comp('HI', zcomp=1., vlim=[-100,100]*u.km/u.s)
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    np.testing.assert_allclose(wvobs_chunks[0][0], 2430.52912749*u.AA)
    np.testing.assert_allclose(wvobs_chunks[0][1], 2432.15114303*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][0], 2050.76022589*u.AA)
    np.testing.assert_allclose(wvobs_chunks[1][1], 2052.12880236*u.AA)
    abscomp, HIlines = mk_comp('HI', zcomp=1., vlim=[-100,100]*u.km/u.s)
    abscomp._abslines[0].analy['wvlim'] = [0,0]*u.AA
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    abscomp._abslines[1].analy['vlim'] = [0,0]*u.AA
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)
    abscomp._abslines[0].attrib['z'] = 0
    abscomp._abslines[0].analy['wvlim'] = [1,0]*u.AA
    wvobs_chunks = ltiu.get_wvobs_chunks(abscomp)


def test_coincident_components():
    abscomp, HIlines = mk_comp('HI', zcomp=2.92939)
    SiIIcomp1,_ = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, zcomp=2.92939)
    SiIIcomp2,_ = mk_comp('SiII',vlim=[-300.,0.]*u.km/u.s, zcomp=2.92939)
    assert ltiu.coincident_components(abscomp, abscomp)  # should overlap
    assert not ltiu.coincident_components(abscomp, SiIIcomp1)  # should not overlap
    assert not ltiu.coincident_components(SiIIcomp2, SiIIcomp1) # should not overlap
    with pytest.raises(ValueError):
        a = ltiu.coincident_components('not_a_component', SiIIcomp1)
    with pytest.raises(ValueError):
        a = ltiu.coincident_components(abscomp, 'not_a_component')


def test_group_coincident_components():
    abscomp, HIlines = mk_comp('HI', zcomp=2.92939)
    SiIIcomp1, _ = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, zcomp=2.92939)
    SiIIcomp2, _ = mk_comp('SiII',vlim=[-300.,0.]*u.km/u.s, zcomp=2.92939)
    # reset names for easy testing
    abscomp.name = 'HI'
    SiIIcomp1.name = 'SiII_1'
    SiIIcomp2.name = 'SiII_2'
    comp_list = [abscomp, abscomp, SiIIcomp1, abscomp, SiIIcomp2, abscomp, SiIIcomp1, SiIIcomp2]
    out = ltiu.group_coincident_components(comp_list)
    assert len(out) == 3  # only three groups
    out_names_0 = [comp.name for comp in out[0]]  # these should be only HI in group 0, and 4 of them
    assert np.sum([n == 'HI' for n in out_names_0]) == 4
    assert len(out[0]) == 4
    out_names = [[],[],[]]
    for ii in range(len(out)):
        for comp in out[ii]:
            out_names[ii] += [comp.name]
    assert out_names == [['HI', 'HI', 'HI', 'HI'], ['SiII_2', 'SiII_2'], ['SiII_1', 'SiII_1']]
    # now a case where are all different
    comp_list = [abscomp, SiIIcomp2 ,SiIIcomp1]
    out = ltiu.group_coincident_components(comp_list)
    for a,b in zip(comp_list, out):
        assert len(b) == 1
        assert a == b[0]
    # check output as dictionary
    out = ltiu.group_coincident_components(comp_list, output_type='dict')
    assert isinstance(out, dict)


def test_synthesize_components():
    #
    SiIIcomp1,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True, skip_synth=True, use_rand=False)
    SiIIcomp1.synthesize_colm(redo_aodm=True)
    #
    SiIIcomp2,_ = mk_comp('SiII',vlim=[50.,300.]*u.km/u.s, add_spec=True, skip_synth=True, use_rand=False)
    SiIIcomp2.synthesize_colm(redo_aodm=True)
    #
    synth_SiII = ltiu.synthesize_components([SiIIcomp1,SiIIcomp2])
    np.testing.assert_allclose(synth_SiII.logN, 13.862454764546792)
    np.testing.assert_allclose(synth_SiII.sig_logN, np.array([0.010146946475971825]*2))
    # Failures
    pytest.raises(IOError, ltiu.synthesize_components, 1)
    pytest.raises(IOError, ltiu.synthesize_components, [1,2])
    # Now SiII*
    SiIIscomp1,_ = mk_comp('SiII*',vlim=[-300.,50.]*u.km/u.s, add_spec=True, skip_synth=True)
    SiIIscomp1.synthesize_colm(redo_aodm=True)
    #
    SiIIscomp2,_ = mk_comp('SiII*',vlim=[50.,300.]*u.km/u.s, add_spec=True, skip_synth=True)
    SiIIscomp2.synthesize_colm(redo_aodm=True)
    synth_SiIIs = ltiu.synthesize_components([SiIIscomp1,SiIIscomp2])
    assert synth_SiIIs.name.count('*') == 1


def test_unique_comps():
    # One list of components
    SiIIcomp,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s)
    HIcomp,_ = mk_comp('HI')
    # Run
    uniq = ltiu.unique_components([HIcomp, SiIIcomp], [SiIIcomp])
    assert np.all(uniq == np.array([True,False]))

def test_add_absline():
    abscomp,_ = mk_comp('HI', zcomp=0.)
    abscomp.add_absline(AbsLine('HI 972'), chk_sep=False, chk_vel=False)
    with pytest.raises(ValueError):
        abscomp.add_absline(AbsLine('HI 949'), vtoler=-10)
    # failed addition
    bad_absline = AbsLine('CIV 1550')
    bad_absline.limits.set([500, 1000]*u.km/u.s)
    bad_absline.attrib['coord'] = SkyCoord(20,20, unit='deg')
    abscomp.add_absline(bad_absline)


def test_fromtodict():
    SiIIcomp1,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    cdict = SiIIcomp1.to_dict()
    #
    assert isinstance(cdict, dict)
    assert cdict['Zion'] == (14, 2)
    # And instantiate
    newcomp = AbsComponent.from_dict(cdict)
    assert isinstance(newcomp, AbsComponent)
    newcomp = AbsComponent.from_dict(cdict, coord=SkyCoord(0,0, unit='deg'))


def test_build_table():
    abscomp,_ = mk_comp('HI')
    # Instantiate
    comp_tbl = abscomp.build_table()
    # Test
    assert isinstance(comp_tbl, Table)
    # empty
    abscomp._abslines = []
    comp_tbl = abscomp.build_table()


def test_synthesize_colm():
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True,
                        add_trans=True, skip_synth=True)
    # Column
    abscomp.synthesize_colm(redo_aodm=True)
    # Test
    np.testing.assert_allclose(abscomp.logN, 13.594445560856554)
    # Reset flags (for testing)
    abscomp2,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True, use_rand=False,
                         add_trans=True, skip_synth=True)
    for iline in abscomp2._abslines:
        if iline.data['name'] == 'SiII 1260':
            iline.attrib['flag_N'] = 2
        elif iline.data['name'] == 'SiII 1808':
            iline.attrib['flag_N'] = 3
        elif iline.data['name'] == 'SiII 1193':
            iline.attrib['flag_N'] = 0
        else:
            iline.attrib['flag_N'] = 1
    abscomp2.synthesize_colm()
    # Test
    np.testing.assert_allclose(abscomp2.logN, 13.3)
    # Another
    abscomp3,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True, use_rand=False, skip_synth=True)
    for iline in abscomp3._abslines:
        if iline.data['name'] == 'SiII 1260':
            iline.attrib['flag_N'] = 3
        elif iline.data['name'] == 'SiII 1808':
            iline.attrib['flag_N'] = 2
        else:
            iline.attrib['flag_N'] = 1
    abscomp3.synthesize_colm()
    # Test
    np.testing.assert_allclose(abscomp3.logN, 13.3)
    # test error
    with pytest.raises(IOError):
        abscomp3.synthesize_colm(overwrite=False)
    with pytest.raises(ValueError):
        abscomp3._abslines[0].attrib['N'] = 0 / u.cm / u.cm
        abscomp3.synthesize_colm(overwrite=True)


def test_build_components_from_lines():
    # Lines
    abscomp,HIlines = mk_comp('HI')
    abscomp,SiIIlines = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # Components
    comps = ltiu.build_components_from_abslines([HIlines[0],HIlines[1],SiIIlines[0],SiIIlines[1]])
    assert len(comps) == 2


def test_abscomp_H2():
    Zion = (-1, -1)  # temporary code for molecules
    Ntuple = (1, 17, 0.5)  # initial guess for Ntuple (needs to be given for adding lines from linelist)
    coord = SkyCoord(0,0, unit='deg')
    z = 0.212
    vlim = [-100., 100.] * u.km/u.s
    comp = AbsComponent(coord, Zion, z, vlim, Ntup=Ntuple)
    comp.add_abslines_from_linelist(llist='H2', init_name="B19-0P(1)", wvlim=[1100, 5000]*u.AA)
    assert len(comp._abslines) == 7


def test_add_abslines_from_linelist():
    comp, HIlines = mk_comp('HI')
    comp._abslines = [] # reset abslines
    comp.add_abslines_from_linelist(llist='HI')
    assert len(comp._abslines) == 30
    comp._abslines = [] # reset
    comp.add_abslines_from_linelist(llist='HI', wvlim=[4100, 5000]*u.AA)
    assert len(comp._abslines) == 1
    # check for no transitions
    comp._abslines = [] # reset
    comp.add_abslines_from_linelist(llist='HI', wvlim=[5000, 5100]*u.AA)
    assert len(comp._abslines) == 0
    # test min_Wr
    comp._abslines = [] # reset
    comp.logN = 13.0
    comp.add_abslines_from_linelist(llist='HI', min_Wr=0.001*u.AA)
    assert len(comp._abslines) == 4
    # test logN not defined
    comp.logN = 0.0
    comp._abslines = [] # reset
    comp.add_abslines_from_linelist(llist='HI', min_Wr=0.001*u.AA)


''' DEPRECATED
def test_iontable_from_components():
    # Lines
    abscomp,HIlines = mk_comp('HI')
    abscomp,SiIIlines = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s, add_spec=True)
    # Components
    comps = ltiu.build_components_from_abslines([HIlines[0],HIlines[1],SiIIlines[0],SiIIlines[1]])
    tbl = ltiu.iontable_from_components(comps)
    assert len(tbl) == 2
'''


def test_complist_from_table_and_table_from_complist():
    # Setup
    tab = Table()
    tab['ion_name'] = ['HI', 'HI', 'CIV', 'SiII', 'OVI']
    tab['Z'] = [1,1,4,14,8]
    tab['ion'] = [1,1,4,2,6]
    tab['z_comp'] = [0.05, 0.0999, 0.1, 0.1001, 0.6]
    tab['RA'] = [100.0]*len(tab) * u.deg
    tab['DEC'] = [-0.8]*len(tab) * u.deg
    tab['vmin'] = [-50.] *len(tab) * u.km / u.s
    tab['vmax'] = [100.] *len(tab) * u.km / u.s
    tab['Ej'] = [0.] *len(tab) / u.cm
    tab['reliability'] = ['a', 'b', 'b', 'none', 'a']
    complist = ltiu.complist_from_table(tab)
    tmp = [-50.,100.]*u.km/u.s
    assert np.all(np.isclose(complist[0].vlim.value, tmp.value))
    tab2 = ltiu.table_from_complist(complist)
    np.testing.assert_allclose(tab['z_comp'], tab2['z_comp'])

    # Without units
    tab = Table()
    tab['ion_name'] = ['HI', 'HI', 'CIV', 'SiII', 'OVI']
    tab['Z'] = [1,1,4,14,8]
    tab['ion'] = [1,1,4,2,6]
    tab['z_comp'] = [0.05, 0.0999, 0.1, 0.1001, 0.6]
    tab['RA'] = [100.0]*len(tab) * u.deg
    tab['DEC'] = [-0.8]*len(tab) * u.deg
    tab['vmin'] = [-50.] *len(tab)
    tab['vmax'] = [100.] *len(tab)
    tab['Ej'] = [0.] *len(tab)
    tab['reliability'] = ['a', 'b', 'b', 'none', 'a']
    complist = ltiu.complist_from_table(tab)
    assert np.all(np.isclose(complist[0].vlim.value, tmp.value))
    assert complist[0].Ej.unit == 1/u.cm

    # test other columns
    tab['logN'] = 13.7
    tab['sig_logN'] = 0.1*np.ones((len(tab),2))
    tab['flag_logN'] = 1
    complist = ltiu.complist_from_table(tab)
    tab2 = ltiu.table_from_complist(complist)
    np.testing.assert_allclose(tab['logN'], tab2['logN'])

    # NHI obj
    SiIIcomp1,_ = mk_comp('SiII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    SiIIIcomp1,_ = mk_comp('SiIII',vlim=[-300.,50.]*u.km/u.s, add_spec=True)
    NHI_obj = DummyObj()
    NHI_obj.NHI = 21.0
    NHI_obj.sig_NHI = np.array([0.1,0.1])
    NHI_obj.flag_NHI = 1
    tab2b = ltiu.table_from_complist([SiIIcomp1, SiIIIcomp1], NHI_obj=NHI_obj)
    assert tab2b['Z'][-1] == 1
    assert tab2b['ion_name'][-1] == 'HI'

    # comment now
    tab['comment'] = ['good', 'good', 'bad', 'bad', 'This is a longer comment with symbols &*^%$']
    complist = ltiu.complist_from_table(tab)
    tab2 = ltiu.table_from_complist(complist)
    comp = complist[-1]
    assert comp.comment == 'This is a longer comment with symbols &*^%$'
    assert tab2['comment'][-1] == comp.comment
    assert tab2['reliability'][-1] == comp.reliability

    # other naming
    tab['name'] = tab['ion_name']
    complist = ltiu.complist_from_table(tab)
    tab2 = ltiu.table_from_complist(complist)
    comp = complist[-1]
    assert comp.name == 'OVI'
    # other attributes
    tab['b'] = [10, 10, 20, 10, 60] * u.km / u.s
    complist = ltiu.complist_from_table(tab)
    tab2 = ltiu.table_from_complist(complist)
    comp = complist[-1]
    assert comp.attrib['b'] == 60*u.km/u.s

    # test errors
    #tab['sig_b'] = [1,2,3,4,5] * u.AA
    #with pytest.raises(IOError):
    #    complist = ltiu.complist_from_table(tab) # bad units for sig_b
    tab = Table()
    tab['ion_name'] = ['HI', 'HI', 'CIV', 'SiII', 'OVI']
    tab['z_comp'] = [0.05, 0.0999, 0.1, 0.1001, 0.6]
    with pytest.raises(IOError):
        complist = ltiu.complist_from_table(tab) # not enough mandatory columns


def test_copy():
    comp, HIlines = mk_comp('HI')
    comp2 = comp.copy()
    comp.logN = 14.6
    assert comp2.logN != comp.logN


def test_add_abslines_from_linelist():
    # Component
    abscomp,_ = mk_comp('SiII', vlim=[-250,80.]*u.km/u.s)
    # Set wavelength (observed) range
    wvrange = [4674. * u.Angstrom, 4690. * u.Angstrom]
    # Add a couple absorption lines
    abscomp.add_abslines_from_linelist(wvlim=wvrange)
    assert len(abscomp._abslines) == 6


""" Utilities for isgm
 Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb
import numpy as np
import warnings
import json

from astropy.table import Table
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord

from linetools import utils as ltu
from linetools.analysis.absline import linear_clm
from linetools.isgm.abssystem import GenericAbsSystem
from linetools.isgm.abscomponent import AbsComponent
from linetools.spectralline import AbsLine
from linetools.lists.linelist import LineList

ckms = const.c.to('km/s').value

def abssys_from_json(filename):
    """
    Parameters
    ----------
    filename

    Returns
    -------
    abs_sys : AbsSystem

    """
    # Load JSON file to determine type
    adict = ltu.loadjson(filename)
    if 'class' in adict.keys():
        if adict['class'] == 'MgIISystem':
            from pyigm.abssys.igmsys import MgIISystem
            abs_sys = MgIISystem.from_dict(adict)
        else:
            warnings.warn("Unknown or uncoded class: {:s}.\nMaking a Generic one".format(adict['class']))
            abs_sys = GenericAbsSystem.from_dict(adict)
    else:
        abs_sys = GenericAbsSystem.from_dict(adict)

    # Return
    return abs_sys


def read_joebvp_to_components(filename, coord, llist=None, specfile=None, chk_vel=False):
    """ Generate a list of AbsComponent objects from a JoeB VP output file

    Parameters
    ----------
    filename : str
      joeB VP filename
    coord : SkyCoord
      QSO sightline
    llist : LineList, optional
      Used to construct AbsLine objects
    specfile : str, optional
    chk_vel : bool, optional
      Demand that the velocities of a given ion all be the same

    Returns
    -------
    comps : list
      list of AbsComponent objects
    """
    # init
    if llist is None:
        llist = LineList('ISM')
    comps = []
    # Read
    vp_data = Table.read(filename, format='ascii')

    # Subset by zsys + trans
    lbls = []
    for izsys, itrans in zip(vp_data['zsys'], vp_data['trans']):
        lbls.append('{:.6f}_{:s}'.format(izsys, itrans))
    lbls = np.array(lbls)
    ulbls = np.unique(lbls)

    # Subset by nflag; Build components
    for lbl in ulbls:
        mt_lines = np.where(lbls == lbl)[0]
        if chk_vel:
            if len(np.unique(vp_data['vel'][mt_lines])) != 1:
                pdb.set_trace()
        z_fit = ltu.z_from_dv(vp_data['vel'][mt_lines[0]]*u.km/u.s, vp_data['zsys'][mt_lines[0]])
        # Loop on abs lines
        alines = []
        for idx in mt_lines:
            zlim = [vp_data['zsys'][idx] +
                    vp_data[vkey][idx] * (1 + vp_data['zsys'][idx]) / ckms
                    for vkey in ['vlim1', 'vlim2']]
            absline = AbsLine(vp_data['restwave'][idx] * u.AA, z=z_fit,
                              zlim=zlim, linelist=llist)
            # Add measurements [JB -- Want to capture anything else??]
            absline.attrib['coord'] = coord
            absline.attrib['flag_N'] = 1
            absline.attrib['logN'] = vp_data['col'][idx]
            absline.attrib['sig_logN'] = vp_data['sigcol'][idx]
            absline.attrib['b'] = vp_data['bval'][idx] * u.km/u.s
            absline.attrib['sig_b'] = vp_data['sigbval'][idx] * u.km/u.s
            absline.attrib['z'] = z_fit
            absline.attrib['sig_z'] = ltu.dz_from_dv(vp_data['sigvel'][idx]*u.km/u.s, vp_data['z_comp'][idx])
            if specfile is None:
                absline.attrib['specfile'] = vp_data['specfile'][idx]
            else:
                absline.attrib['specfile'] = specfile
            # Fill N, sig_N
            _, _, = linear_clm(absline.attrib)
            alines.append(absline)

        # AbsComponent
        stars = '*' * alines[0].ion_name.count('*')
        if 'comment' in vp_data.keys():
            comment = vp_data['comment'][mt_lines[0]]
        else:
            comment = ''
        if 'rely' in vp_data.keys():
            reliability = vp_data['rely'][mt_lines[0]]
        else:
            reliability = 'none'
        abscomp = AbsComponent.from_abslines(alines, stars=stars, comment=comment, reliability=reliability)

        # Add measurements [JB -- Want to capture anything else??]
        abscomp.attrib = alines[0].attrib.copy()
        # Remove undesired keys
        for key in ['EW', 'sig_EW', 'flag_EW', 'N', 'sig_N']:
            abscomp.attrib.pop(key)
        # And more required
        for key in ['flag_N', 'logN', 'sig_logN']:
            setattr(abscomp, key, abscomp.attrib[key])
        # Errors must be in first line!
        assert abscomp.sig_logN > 0., "AbsComponent has sig_logN=0 {}".format(abscomp)

        comps.append(abscomp)
    # Finish
    return comps


def write_joebvp_from_components(comp_list, specfile, outfile):
    """ From a given component list, it produces an
    input file for JOEBVP (Voigt profile fitter).

    Parameters
    ----------
    comp_list : list of AbsComponent
        Input list of components to group
    specfile : str
        Name of the spectrum file associated to the components
        in comp_list
    outfile : str
        Name of the output file

    """
    # Open new file to write out
    f = open(outfile, 'w')

    # Print header
    s = 'specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|z_comp|trans|rely|comment\n'
    f.write(s)

    # Components
    for ii, comp in enumerate(comp_list):
        flags = (ii+2,ii+2,ii+2)
        try:
            b_val = comp.attrib['b']
        except KeyError:
            b_val = 10*u.km/u.s
        s = comp.repr_joebvp(specfile, flags=flags, b_default=b_val)  # still, b values from abslines take precedence if they exist
        f.write(s)
    f.close()


def read_igmg_to_components(igmguesses_json):
    """Creates a list of AbsComponenbts from a igmguesses_json file

    Parameters
    ----------
    igmguesses_json : str
        Name of the json file generated by IGMGuesses

    Returns
    -------
    comps : list
        A list of AbsComponents

    """
    # Read the JSON file
    with open(igmguesses_json) as data_file:
        igmg_dict = json.load(data_file)
    # Components
    comp_list = []
    icoord = SkyCoord(ra=igmg_dict['meta']['RA'], dec=igmg_dict['meta']['DEC'], unit='deg')
    for ii, key in enumerate(igmg_dict['cmps'].keys()):
        # import pdb;pdb.set_trace()
        comp_dict = igmg_dict['cmps'][key]
        comp_dict['flag_N'] = 1
        comp_dict['logN'] = comp_dict['Nfit']
        comp_dict['sig_logN'] = -1
        try:
            comp = AbsComponent.from_dict(comp_dict, coord=icoord, chk_sep=False, chk_data=False, chk_vel=False, linelist="ISM")
        except:  # is this a H2 line?
            comp = AbsComponent.from_dict(comp_dict, coord=icoord, chk_sep=False, chk_data=False, chk_vel=False, linelist='H2')
        comp_list += [comp]
    return comp_list


def write_igmg_from_components(comp_list, specfile, fwhm, outfile='IGM_model.json'):
        """ Write to a JSON file of the IGMGuesses format from a list of AbsComponents.
        All coordinates must be equal!

        Parameters
        ----------
        comp_list : list
            list of AbsComponents
        specfile : str
            Name of spectrum file associated to these components
        fwhm : int
            FWHM of the spectrum
        outfile : str, optional
            Name of the output json file

        """
        import datetime
        # coordinates and meta
        coord_ref = comp_list[0].coord
        RA = coord_ref.ra.to('deg').value
        DEC = coord_ref.dec.to('deg').value
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        jname = ltu.name_from_coord(coord_ref, precision=(2, 1))
        # Create dict of the components
        out_dict = dict(cmps={},
                        spec_file=specfile,
                        fwhm=fwhm, bad_pixels=[],
                        meta={'RA': RA, 'DEC': DEC, 'ALTNAME':'unknown',
                              'zem': 0., 'Creator': 'linetools from AbsComponent list',
                              'Instrument': 'unknown', 'Date': date, 'JNAME': jname})

        for comp in comp_list:
            key = comp.name
            out_dict['cmps'][key] = comp.to_dict()
            # import pdb; pdb.set_trace()
            # check coordinate
            if comp.coord != coord_ref:
                raise ValueError("All AbsComponent objects must have the same coordinates!")
            out_dict['cmps'][key]['zcomp'] = comp.zcomp
            out_dict['cmps'][key]['zfit'] = comp.zcomp
            out_dict['cmps'][key]['Nfit'] = comp.logN
            out_dict['cmps'][key]['bfit'] = comp.attrib['b']
            out_dict['cmps'][key]['wrest'] = comp._abslines[0].wrest.value
            out_dict['cmps'][key]['vlim'] = list(comp.vlim.value)
            out_dict['cmps'][key]['reliability'] = str(comp.reliability)
            out_dict['cmps'][key]['comment'] = str(comp.comment)
            # out_dict['cmps'][key]['mask_abslines'] = comp.mask_abslines

        # JSONify
        gd_dict = ltu.jsonify(out_dict)

        # Write file
        # with io.open(outfile, 'w', encoding='utf-8') as f:
        f = open(outfile, 'w')
        f.write(unicode(json.dumps(gd_dict, sort_keys=True, indent=4, separators=(',', ': '))))
        print('Wrote: {:s}'.format(outfile))
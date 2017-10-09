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

from astropy.table import Table
from astropy import constants as const
from astropy import units as u

from linetools import utils as ltu
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


def read_joebvp(filename, coord, llist=None):
    """
    Parameters
    ----------
    filename : str
      joeB VP filename
    coord : SkyCoord
    llist : LineList, optional

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
    nflags, nidx = np.unique(vp_data['nflag'], return_index=True)

    # Subset by nflag; Build components
    for nflag in nflags:
        mt_lines = np.where(vp_data['nflag'] == nflag)[0]
        if len(np.unique(vp_data['vel'][mt_lines])) != 1:
            pdb.set_trace()
        z_fit = vp_data['zsys'][mt_lines[0]] + vp_data['vel'][mt_lines[0]] * (1 + vp_data['zsys'][mt_lines[0]]) / ckms
        # Loop on abs lines
        alines = []
        for idx in mt_lines:
            zlim = [vp_data['zsys'][idx] +
                    vp_data[vkey][idx] * (1 + vp_data['zsys'][idx]) / ckms
                    for vkey in ['vlim1', 'vlim2']]
            absline = AbsLine(vp_data['restwave'][idx] * u.AA, z=z_fit, zlim=zlim, linelist=llist)
            # Add measurements [JB -- Want to capture anything else??]
            absline.attrib['coord'] = coord
            absline.attrib['flag_N'] = 1
            absline.attrib['logN'] = vp_data['col'][idx]
            absline.attrib['sig_logN'] = vp_data['sigcol'][idx]
            absline.attrib['b'] = vp_data['bval'][idx]
            absline.attrib['sig_b'] = vp_data['sigbval'][idx]
            absline.attrib['z'] = z_fit
            absline.attrib['sig_z'] = (1 + vp_data['z_comp'][idx]) * vp_data['sigvel'][idx] / ckms
            absline.attrib['specfile'] = vp_data['specfile'][idx]
            alines.append(absline)
        # Component
        stars = '*' * alines[0].ion_name.count('*')
        abscomp = AbsComponent.from_abslines(alines, stars=stars)
        # Add measurements [JB -- Want to capture anything else??]
        abscomp.attrib = alines[0].attrib.copy()
        # Remove undesired keys
        for key in ['EW', 'sig_EW', 'flag_EW', 'N', 'sig_N']:
            abscomp.attrib.pop(key)
        # And more
        for key in ['flag_N', 'logN', 'sig_logN']:
            setattr(abscomp, key, abscomp.attrib[key])
        # Errors must be in first line!
        assert abscomp.sig_logN > 0.
        comps.append(abscomp)
    # Finish
    return comps
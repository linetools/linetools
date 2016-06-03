""" Module for utilities related to spectral line or lines(s)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np

from astropy.table import Table
from astropy.units import Quantity

from .spectralline import AbsLine

def parse_speclines(abslines, key, mk_array=False):
    """ Generate a list or array of items from a list of AbsLines

    Parameters
    ----------
    abslines : list of AbsLine objects
    key : str

    Returns
    -------
    items : list

    """
    out_list = []
    # Ugly loop
    for iline in abslines:
        try:
            out_list.append(getattr(iline,key))
        except AttributeError:
            try:
                out_list.append(iline.attrib[key])
            except KeyError:
                try:
                    out_list.append(iline.analy[key])
                except KeyError:
                    try:
                        out_list.append(iline.data[key])
                    except KeyError:
                        out_list.append(None)
    # Return
    if mk_array:
        if isinstance(out_list[0], Quantity):
            return Quantity(out_list)
        else:
            return np.array(out_list)
    else:
        return out_list


def transtable_from_speclines(abslines, add_keys=None):
    """Generate a Table summarizing the transitions from a list of AbsLines
    Parameters
    ----------
    speclines : list of SpectralLine objects

    Returns
    -------
    tbl : Table

    """
    keys = ['wrest','name','Z', 'ion', 'Ej', 'z', 'EW', 'sig_EW']
    if isinstance(abslines[0], AbsLine):
        keys += ['flag_N', 'logN', 'sig_logN']
    if add_keys is not None:
        keys += add_keys

    # Get started
    tbl = Table()

    # Loop to my loop
    for key in keys:
        tbl[key] = parse_speclines(abslines, key, mk_array=True)

    # Sort
    tbl.sort(['Z','ion','Ej','wrest'])
    # Return
    return tbl

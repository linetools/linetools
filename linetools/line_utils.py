""" Module for utilities related to spectral line or lines(s)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np

from astropy.table import QTable, Column
from astropy import units as u
from astropy.units import Quantity

def parse_abslines(abslines, key, mk_array=False):
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


def transtable_from_abslines(abslines, add_keys=None):
    """Generate a QTable summarizing the transitions from a list of AbsLines
    Parameters
    ----------
    abslines : list of AbsLine objects

    Returns
    -------
    tbl : QTable

    """
    keys = ['wrest','name','z','flag_EW', 'EW', 'flag_N', 'logN']
    if add_keys is not None:
        keys += add_keys

    # Get started
    tbl = QTable()

    # Loop to my loop
    for key in keys:
        tbl[key] = parse_abslines(abslines, key, mk_array=True)

    # Return
    return tbl

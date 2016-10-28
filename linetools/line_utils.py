""" Module for utilities related to spectral line or lines(s)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np

from astropy.table import Table
from astropy.units import Quantity

from linetools.spectralline import AbsLine, EmLine

def parse_speclines(speclines, key, mk_array=False):
    """ Generate a list or array of items from a list of SpectralLines
    for a given key

    Parameters
    ----------
    abslines : list of AbsLine objects
    key : str
      Property of interest, e.g. 'wrest', 'EW', 'f'
    mk_array : bool, optional
      Return an array (or Quantity array) instead of a list

    Returns
    -------
    items : list or array

    """
    out_list = []
    # Ugly loop
    for iline in speclines:
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


def transtable_from_speclines(speclines, add_keys=None):
    """Generate a Table summarizing the transitions from a list of SpectralLines

    Parameters
    ----------
    speclines : list of SpectralLine objects
    add_keys : list, optional
      Additional keys to include in Table

    Returns
    -------
    tbl : Table

    """
    keys = ['wrest','name','Z', 'ion', 'Ej', 'z', 'EW', 'sig_EW']
    if speclines[0].ltype == 'Abs':
        keys += ['flag_N', 'logN', 'sig_logN']
    if speclines[0].ltype == 'Em':
        keys += ['flag_flux', 'flux', 'sig_flux']
    if add_keys is not None:
        keys += add_keys

    # Get started
    tbl = Table()

    # Loop to my loop
    for key in keys:
        tbl[key] = parse_speclines(speclines, key, mk_array=True)

    # Sort
    tbl.sort(['Z','ion','Ej','wrest'])
    # Return
    return tbl

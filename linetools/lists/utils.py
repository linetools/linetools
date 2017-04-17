""" Utilities related to line lists
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

from astropy.table import QTable, Table
from astropy import units as u

def from_dict_to_table(a):
    """Converts dictionary `a` to its Table version.

    Parameters
    ----------
    a : dict
        The input dictionary to be converted to a
        Table. The resulting Table will have the
        same keys as self._data

    Returns
    -------
    A Table of 1 Row, with filled with the data from
    the input dictionary.

    """
    # Check
    if isinstance(a, dict):
        pass
    else:
        raise ValueError('Input has to be a dictionary.')
    #
    tdict = {}
    for key in a.keys():
        if isinstance(a[key], u.Quantity):
            tdict[key] = [a[key].value] * a[key].unit
        else:
            tdict[key] = [a[key]]
    # Table
    tab = Table(tdict)
    # Return
    return tab


def from_table_to_dict(tab):
    """Converts Table `tab` to its dictionary version.
    Adds units when appropriate
    An error is raised if len(tab) > 1.

    Parameters
    ----------
    tab : QTable or Table
        The table to be converted to a dictionary.
        It has to be of length == 1!

    Returns
    -------
    A dictionary with the same keys and values of the
    input table.

    """

    if not isinstance(tab, (QTable, Table)):
        raise ValueError('Input has to be QTable or Table.')
    elif len(tab) != 1:
        raise ValueError('Input has to be of len(tab) == 1.')

    a_dict = dict()
    for key in tab.keys():
        # Units
        if tab[key].unit is not None:
            a_dict[key] = tab[0][key] * tab[key].unit
        else:
            a_dict[key] = tab[0][key]
    return a_dict
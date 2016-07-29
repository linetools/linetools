""" Utilities for working with ionized atoms.
"""
#;+
#; NAME:
#; ionization
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for ionization of atoms
#;   03-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

from linetools.abund.elements import ELEMENTS
from linetools.abund import roman

########################## ##########################
########################## ##########################
def ion_name(ion, flg=0, nspace=None):
    """ Convert ion tuple into a string

    Parameters
    ----------
    ion : tuple or dict
      Either a tuple of integers (Z, ion) or a dict with tags of `Z`
      and `ion`. e.g. (6, 4) would return 'CIV'.

    flg : int, optional (0)
        * 0: Roman numeral (e.g. CIV)
        * 1: Latex with ion notation (e.g C^+)

    nspace : int, optional (0)
      Number of spaces to insert.

    Returns
    -------
    name : str
      e.g. SiII, {\\rm Si}^{+}

    """
    if isinstance(ion,tuple):
        elm = ELEMENTS[ion[0]]
        str_elm = elm.symbol
    else:
        return ion_name( (ion['Z'], ion['ion']) )

    # Ion state
    if flg == 0: # Roman
        if nspace is None:
            nspace = 0
        str_ion = roman.toRoman(ion[1])
        spc = ' ' * nspace
        outp = str_elm + spc + str_ion
    elif flg == 1: # LaTeX
        if ion[1] == 0:
            raise ValueError('ionization.ion_name: Not ready for this input yet.')
        elif ion[1] == 1:
            str_ion = '^0'
        elif ion[1] == 2:
            str_ion = '^{+}'
        elif ion[1] == 3:
            str_ion = '^{++}'
        else:
            str_ion = '^{+' + str(ion[1] - 1) + '}'
        outp = '{\\rm ' + str_elm + '}' + str_ion
    else:
        raise ValueError('ionization.ion_name: Not ready for this flg.')

    return outp


########################## ##########################
########################## ##########################
def name_ion(ion):
    """ Convert string into ion tuple

    Parameters
    ----------
    ion : str
      Name of the ion, e.g. 'SiII' or 'Si II'

    Returns
    -------
    ion_tup : tuple
      Z, ion -- e.g. (14,2)
    """
    if isinstance(ion,basestring):
        pass
    else:
        raise ValueError('ionization.name_ion: Not ready for this input yet.')

    ion = ion.strip('*') # e.g. CII*

    if ion[1] in ['I','V', 'X', ' ']:
        iion = 1
    else:
        iion = 2

    # Element
    elm = ion[0:iion]
    if elm == 'D': # Deuterium
        Z = 1
    else:
        Z = ELEMENTS[ion[0:iion]].number

    # Ion
    ion_state = roman.fromRoman(ion[iion:].strip())

    return Z, ion_state

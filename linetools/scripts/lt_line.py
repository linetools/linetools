#!/usr/bin/env python

"""
Plot the data for a line (or lines)
Also print the line data
  Examples:
  lt_line HI
  lt_line HI 1215
  lt_line 1215
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:
    ustr = unicode
except NameError:
    ustr = str

def main(args=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Print spectral line data of a line or lines.')
    parser.add_argument("inp1", nargs='?', default=None, help="Ion, transition name, or Rest wavelength (e.g. HI, HI1215 or 1215)")
    parser.add_argument("-all", default=False, action='store_true', help="Print all lines")
    parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")
    parser.add_argument("-toler", default=1., action='store_true', help="Matching tolerance (in Ang) on input wavelength")
    #parser.add_argument("wave", type=str, default=None, help="Rest wavelength in Angstroms")

    pargs = parser.parse_args()
    if pargs.inp1 is None and pargs.all is False:
        print("No option selected.  Use -h for Help")
        return
    # Setup
    from astropy import units as u
    from linetools.lists.linelist import LineList
    import numpy as np
    cols = ['name', 'wrest', 'f', 'A']
    # LineList
    llist = LineList(pargs.llist)
    # All?
    if pargs.all:
        llist._data[cols].pprint(99999)
        return
    # Grab line(s)
    if ustr(pargs.inp1[0]).isdecimal():
        wrest = float(pargs.inp1)*u.AA
        mtch = np.abs(wrest-llist.wrest) < pargs.toler*u.AA
        llist._data[cols][mtch].pprint(99999)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

"""
Print the solar abundance data for an element or all elements
  Examples:
  lt_solabnd Fe
  lt_solabnd -a
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Print Solar abundance data for an element or all elements.')
    parser.add_argument("inp", nargs='?', default=None, help="Elm (e.g. H, Fe)")
    parser.add_argument("-a", "--all", default=False, action='store_true', help="Print all values")
    parser.add_argument("--sortZ", default=False, action='store_true', help="Sort on Atomic Number")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    pargs = parser(options=args)
    if pargs.inp is None and pargs.all is False:
        print("No option selected.  Use -h for Help")
        return
    # Setup
    from astropy import units as u
    from astropy.table import Column
    from linetools.abund.solar import SolarAbund
    from linetools.abund import ions as ltai
    import numpy as np

    print('-----------------------')
    sol = SolarAbund()
    print('-----------------------')
    # All?
    if pargs.all:
        if pargs.sortZ:
            sol._data.sort('Z')
        else:
            sol._data.sort('Elm')
        sol._data.pprint(max_lines=1000)
    else:  # Input better have been set
        imt = np.where(sol._data['Elm'] == pargs.inp)[0]
        if len(imt) == 0:
            print("Bad Element name.  Try again")
        else:
            print(sol._data[imt])

if __name__ == '__main__':
    main()

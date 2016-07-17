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

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Print spectral line data of a line or lines.')
    parser.add_argument("inp", nargs='?', default=None, help="Ion, transition name, or Rest wavelength (e.g. HI, HI1215 or 1215)")
    parser.add_argument("-a", "--all", default=False, action='store_true', help="Print all lines")
    parser.add_argument("--llist", default='ISM', type=str, help="Name of LineList:  ISM, HI, H2, CO, etc.")
    parser.add_argument("-t", "--toler", default=1., type=float, help="Matching tolerance (in Ang) on input wavelength")
    parser.add_argument("-z", "--redshift", default=0., type=float, help="Matching tolerance (in Ang) on input wavelength")
    #parser.add_argument("wave", type=str, default=None, help="Rest wavelength in Angstroms")

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
    from linetools.lists.linelist import LineList
    from linetools.abund import ions as ltai
    import numpy as np
    if pargs.llist in ['CO']:
        cols = ['mol', 'label', 'wrest', 'f']
    elif pargs.llist in ['H2']:
        cols = ['mol', 'name', 'wrest', 'f']
    else:
        cols = ['name', 'wrest', 'f', 'A']
    # LineList
    llist = LineList(pargs.llist)
    # Redshift?
    if float(pargs.redshift) != 0.:
        z = llist._data['wrest']*(1+float(pargs.redshift))
        llist._data.add_column(Column(z, name='z'))
        cols += ['z']
    # All?
    if pargs.all:
        try:
            llist._data[cols].pprint(99999)
        except ValueError:
            pdb.set_trace()
        return
    # Grab line(s)
    if ustr(pargs.inp[0]).isdecimal():  # Input rest wavelength
        wrest = float(pargs.inp)*u.AA
        mtch = np.abs(wrest-llist.wrest) < pargs.toler*u.AA
        llist._data[cols][mtch].pprint(99999)
    else:  # Either ion or transition
        istrans = False
        for jj,char in enumerate(pargs.inp):
            if char.isdigit():
                istrans = True
                i0 = jj
                break
        if istrans:
            trans = pargs.inp[0:i0]+' '+pargs.inp[i0:]
            tdict = llist[trans]
            for key,value in tdict.items():
                if key in cols:
                    print('{:s}: {}'.format(key,value))
        else:  # Ion
            Zion = ltai.name_ion(pargs.inp)
            mtion = (llist.Z == Zion[0]) & (llist.ion == Zion[1])
            llist._data[cols][mtion].pprint(99999)


if __name__ == '__main__':
    main()

#!/usr/bin/env python

"""
Ingest input coordinates and print them to the screen
  Examples:
  lt_radec 152.25900,7.22885
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
        description='Print coordinates in several formats from input one.')
    parser.add_argument("inp", nargs='?', default=None, help="RA,DEC (e.g. 152.25900,7.22885), JXX (e.g. J100902.16+071343.8)")
    parser.add_argument("--epoch", default=2000., type=float, help="Epoch [Not functional]")

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
    from linetools import utils as ltu
    from .utils import coord_arg_to_coord
    from astropy import units as u

    # RA, DEC
    icoord = coord_arg_to_coord(pargs.inp)
    coord = ltu.radec_to_coord(icoord)

    # Time to print
    print('J{:s}{:s}'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True),
                             coord.dec.to_string(sep='',pad=True,alwayssign=True)))
    print('   RA={:f} deg, DEC={:f} deg'.format(coord.ra.deg, coord.dec.deg))
    print('   radec = ({:f},{:f}) deg'.format(coord.ra.deg, coord.dec.deg))
    print('SDSS finding chart: https://skyserver.sdss.org/dr12/en/tools/chart/navi.aspx?ra={:f}&dec={:f}&opt='.format(coord.ra.deg, coord.dec.deg))

if __name__ == '__main__':
    main()

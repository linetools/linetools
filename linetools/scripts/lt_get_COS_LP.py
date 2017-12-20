#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division, unicode_literals)


try:
    ustr = unicode
except NameError:
    ustr = str

"""Obtain the COS life-time position given a date.
"""
def parser(options=None):
    import argparse
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Obtain the HST/COS life-time position (LP) given a date.')
    parser.add_argument("date", type=str, default=None, help='Enter a date, e.g. "2017-10-09 00:00:00", "2017-10-09", etc.')
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args=None):
    from linetools.spectra import utils as lsu

    pargs = parser(options=args)
    date = pargs.date

    # return
    lp = lsu.get_COS_LP_from_date(date)
    print(lp)

if __name__ == '__main__':
    main()
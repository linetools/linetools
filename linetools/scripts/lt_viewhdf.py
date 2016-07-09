#!/usr/bin/env python

"""
View the contents of an HDF5 file
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import h5py

import pdb

def show(obj, level, full=False):
    """ Simple recursive show def

    Parameters
    ----------
    obj

    Returns
    -------
    """
    tlevel = level + 1
    pad = ' ' + '  '*level
    for key in obj.keys():
        if isinstance(obj[key], h5py._hl.group.Group):
            grps = '{:s}Group: {:s}'.format(pad, key)
            print('{:s}'.format(pad)+'-'*(len(grps)-len(pad)))
            print(grps)
            print('{:s}'.format(pad)+'-'*(len(grps)-len(pad)))
            show(obj[key], tlevel, full=full)
        elif isinstance(obj[key], h5py._hl.dataset.Dataset):
            ss = '{:s}{:s}:  {:s}'.format(pad, key, str(obj[key]))
            print('{:s}'.format(pad)+'='*(len(ss)-len(pad)))
            print(ss)
            print('{:s}'.format(pad)+'='*(len(ss)-len(pad)))
            if full:
                nms = obj[key].dtype.descr
                for nm in nms:
                    print('{:s}    {:s}'.format(pad, nm))
            else:
                nms = obj[key].dtype.names
                print('{:s}    {}'.format(pad, nms))
        else:
            raise ValueError("Uh oh")
    #
    return


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='View contents of an HDF5 file')
    parser.add_argument("file", type=str, help="Filename")
    parser.add_argument("-f", "--full", default=False, help="Full output?", action='store_true')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    pargs = parser(options=args)
    # Setup
    # Load
    print("Opening {:s}".format(pargs.file))
    hf5 = h5py.File(pargs.file, 'r')
    # Show
    show(hf5, 0, full=pargs.full)


if __name__ == '__main__':
    main()

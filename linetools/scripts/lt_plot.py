#!/usr/bin/env python
""" Plot a spectrum
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

def plotspec(args):
    """Plot spectrum files
    """
    from linetools.spectra.io import readspec
    import warnings
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    warnings.simplefilter('ignore', mpl.mplDeprecation)

    plt.rcParams['axes.formatter.useoffset'] = False  # avoid scientific notation in axes tick labels

    spec_cache = {}

    fig = plt.figure(figsize=(10,5))
    fig.subplots_adjust(left=0.07, right=0.95, bottom=0.11)
    ax = fig.add_subplot(111)
    i = 0
    quit = False
    print("#### Use left and right arrow keys to navigate, 'Q' to quit ####")

    while quit is False:
        filename = args.filenames[i]
        if filename not in spec_cache:
            spec_cache[filename] = readspec(filename)
        sp = spec_cache[filename]
        ax.cla()
        sp.plot(show=False)
        ax.set_xlabel(str(sp.wavelength.unit))
        ax.set_title(filename)
        if args.redshift is not None:
            from linetools.lists.linelist import LineList
            ll = LineList('Strong')
            #import pdb ;pdb.set_trace()
            wlines = ll._data['wrest'] * (1 + args.redshift)
            y0, y1 = ax.get_ylim()
            ax.vlines(wlines.to(sp.wavelength.unit).value, y0, y1,
                      linestyle='dotted')

        while True:
            plt.waitforbuttonpress()
            if sp._plotter.last_keypress == 'right':
                i += 1
                i = min(i, len(args.filenames) - 1)
                # Note this only breaks out of the inner while loop
                break
            elif sp._plotter.last_keypress == 'left':
                i -= 1
                i = max(i, 0)
                break
            elif sp._plotter.last_keypress == 'Q':
                quit = True
                break


def main(args=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Plot spectrum files')
    parser.add_argument("filenames", nargs='+',
                        help="Spectrum filenames")
    parser.add_argument("-z", "--redshift", type=float, default=None,
                        help="Redshift for line plotting")

    args = parser.parse_args()
    plotspec(args)

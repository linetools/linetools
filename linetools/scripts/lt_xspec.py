#!/usr/bin/env python

"""
Plot a spectrum with an interactive QT GUI
"""
import pdb
import sys


# Script to run XSpec from the command line or ipython
def main(*args, **kwargs):
    """ Runs the XSpecGui on an input file
    """
    import argparse

    parser = argparse.ArgumentParser(description='Parser for lt_xspec v1.2; \n Note: Extra arguments are passed to read_spec (e.g. --flux_tag=FX)')
    parser.add_argument("file", type=str, help="Spectral file; specify extension by appending #exten#")
    parser.add_argument("-guessfile", "--guessfile", type=str, help="Igmguesses file, see https://github.com/pyigm/pyigm/blob/master/docs/igmguesses.rst ")
    parser.add_argument("-z", "--zsys", type=float, help="System Redshift")
    parser.add_argument("--norm", help="Show spectrum continuum normalized (if continuum is provided)",
                        action="store_true")
    parser.add_argument("--air", default=False, help="Convert input spectrum wavelengths from air to vacuum", action="store_true")
    parser.add_argument("--exten", type=int, help="FITS extension")
    parser.add_argument("--splice", type=str, help="Splice with the input file; extension convention applies")
    parser.add_argument("--scale", type=float, help="Scale factor for GUI size [1. is default]")
    #parser.add_argument("--wave_tag", type=str, help="Tag for wave in Table")
    #parser.add_argument("--flux_tag", type=str, help="Tag for flux in Table")
    #parser.add_argument("--sig_tag", type=str, help="Tag for sig in Table")
    #parser.add_argument("--var_tag", type=str, help="Tag for var in Table")
    #parser.add_argument("--ivar_tag", type=str, help="Tag for ivar in Table")

    #pargs = parser.parse_args()
    pargs, unknown = parser.parse_known_args()

    from PyQt5.QtWidgets import QApplication
    from linetools.guis.xspecgui import XSpecGui

    # Normalized?
    if pargs.norm is True:
        norm = True
    else:
        norm = False

    # Extension
    file = pargs.file
    if pargs.file[-1] == '#':
        prs = pargs.file.split('#')
        exten = int(prs[1])
        file = prs[0]
    else:
        exten = (pargs.exten if hasattr(pargs, 'exten') else 0)

    # zsys
    zsys = (pargs.zsys if hasattr(pargs, 'zsys') else None)

    # guesses
    guessfile = (pargs.guessfile if hasattr(pargs, 'guessfile') else None)


    # Splice?
    if pargs.splice is not None:
        pdb.set_trace()
        if pargs.splice[-1] == '#':
            prs = pargs.splice.split('#')
            exten = [exten, int(prs[1])]
            file = [file, prs[0]]
        else:
            exten = [exten, None]
            file = [file, pargs.splice]


    # Read spec keywords
    rsp_kwargs = {}
    for arg in unknown:
        spl = arg.split('=')
        rsp_kwargs[spl[0][2:]] = spl[1]

    # GUI
    app = QApplication(sys.argv)

    # Scale
    if pargs.scale is None:
        # Screen dimensions
        width = app.desktop().screenGeometry().width()
        scale = 2. * (width/3200.)
    else:
        scale = pargs.scale
    #
    gui = XSpecGui(file, guessfile=guessfile, zsys=zsys, norm=norm, exten=exten,
                   rsp_kwargs=rsp_kwargs, air=pargs.air,
                   screen_scale=scale)
    gui.show()
    app.exec_()

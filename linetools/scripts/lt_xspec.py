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

    parser = argparse.ArgumentParser(description='Parse for XSpec')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--zsys", type=float, help="System Redshift")
    parser.add_argument("--norm", help="Show spectrum continuum normalized (if one exists)",
                        action="store_true")
    parser.add_argument("--exten", type=int, help="FITS extension")

    pargs = parser.parse_args()


    from PyQt4 import QtGui
    from linetools.guis.xspecgui import XSpecGui

    # Normalized?
    if pargs.norm is True:
        norm = True
    else:
        norm = False

    # Extension
    exten = (pargs.exten if hasattr(pargs, 'exten') else 0)

    # Second spectral file?
    zsys = (pargs.zsys if hasattr(pargs, 'zsys') else None)


    app = QtGui.QApplication(sys.argv)
    gui = XSpecGui(pargs.file, zsys=zsys, norm=norm, exten=exten)
    gui.show()
    app.exec_()

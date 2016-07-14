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
    parser.add_argument("--norm", help="Show spectrum continuum normalized (if continuum is provided)",
                        action="store_true")
    parser.add_argument("--air", default=False, help="Convert input spectrum wavelengths from air to vacuum", action="store_true")
    parser.add_argument("--exten", type=int, help="FITS extension")
    parser.add_argument("--wave_tag", type=str, help="Tag for wave in Table")
    parser.add_argument("--flux_tag", type=str, help="Tag for flux in Table")
    parser.add_argument("--sig_tag", type=str, help="Tag for sig in Table")
    parser.add_argument("--var_tag", type=str, help="Tag for var in Table")
    parser.add_argument("--ivar_tag", type=str, help="Tag for ivar in Table")

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

    # Read spec keywords
    rsp_kwargs = {}
    if pargs.wave_tag is not None:
        rsp_kwargs['wave_tag'] = pargs.wave_tag
    if pargs.flux_tag is not None:
        rsp_kwargs['flux_tag'] = pargs.flux_tag
    if pargs.sig_tag is not None:
        rsp_kwargs['sig_tag'] = pargs.sig_tag
    if pargs.var_tag is not None:
        rsp_kwargs['var_tag'] = pargs.var_tag
    if pargs.ivar_tag is not None:
        rsp_kwargs['ivar_tag'] = pargs.ivar_tag

    app = QtGui.QApplication(sys.argv)

    gui = XSpecGui(pargs.file, zsys=zsys, norm=norm, exten=exten,
                   rsp_kwargs=rsp_kwargs, air=pargs.air)
    gui.show()
    app.exec_()

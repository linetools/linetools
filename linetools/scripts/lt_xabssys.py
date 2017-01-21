#!/usr/bin/env python

"""
Show a VelocityPlot to interactively modify lines in an AbsSystem
"""
import pdb
import sys


# Script to run XAbsSysGui from the command line or ipython
def main(*args, **kwargs):
    """ Runs the XAbsSysGui on input files
    """
    import argparse

    parser = argparse.ArgumentParser(description='Parse for XAbsSys')
    parser.add_argument("spec_file", type=str, help="Spectral file")
    parser.add_argument("abssys_file", type=str, help="AbsSys file (JSON)")
    parser.add_argument("-outfile", type=str, help="Output filename")
    parser.add_argument("-llist", type=str, help="Name of LineList")
    #parser.add_argument("-exten", type=int, help="FITS extension")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized",
                        action="store_true")

    pargs = parser.parse_args()

    from PyQt4 import QtGui
    from linetools.guis.xabssysgui import XAbsSysGui

    # Normalized?
    norm = True
    if pargs.un_norm:
        norm = False

    # Extension
    #exten = (pargs.exten if hasattr(pargs, 'exten') else 0)

    # Read spec keywords
    rsp_kwargs = {}

    # Line list
    if pargs.llist is not None:
        from linetools.lists.linelist import LineList
        llist = LineList(pargs.llist)
    else:
        llist = None

    # Read AbsSystem
    from linetools.isgm.abssystem import GenericAbsSystem
    abs_sys = GenericAbsSystem.from_json(pargs.abssys_file)#, chk_vel=False)

    app = QtGui.QApplication(sys.argv)

    gui = XAbsSysGui(pargs.spec_file, abs_sys, norm=norm, llist=llist,
                     outfil=pargs.outfile)
    gui.show()
    app.exec_()

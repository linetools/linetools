#!/usr/bin/env python

"""
Show a VelocityPlot to interactively modify lines in an AbsSystem
"""
import sys
import pdb
import warnings


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
    parser.add_argument("--specdb", help="Spectral file is a SPECDB database", action="store_true")
    parser.add_argument("--group", type=str, help="SPECDB group name")
    parser.add_argument("--un_norm", help="Spectrum is NOT normalized", action="store_true")
    parser.add_argument("--chk_z",  help="Check the z limits of your components? [default=False]",
                        action="store_true")

    pargs = parser.parse_args()

    from PyQt5.QtWidgets import QApplication
    from linetools.guis.xabssysgui import XAbsSysGui
    from linetools.isgm.io import abssys_from_json
    from IPython import embed

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
    if not pargs.chk_z:
        warnings.warn("Not checking your system's velocity limits.  This is the Default but be so warned.")
    abs_sys = GenericAbsSystem.from_json(pargs.abssys_file, chk_z=pargs.chk_z)
    if len(abs_sys.list_of_abslines()) == 0:
        warnings.warn("No absorption lines given.  I hope you intended that to be the case!")

    app = QApplication(sys.argv)

    # Load spectrum using specdb?
    if pargs.specdb:
        # Instantiate
        from specdb.specdb import SpecDB
        from specdb import group_utils
        sdb = SpecDB(db_file=pargs.spec_file)
        # Grab spectrum
        if pargs.group is not None:
            groups = [pargs.group]
        else:
            groups = None
        spec, meta = sdb.spectra_from_coord(abs_sys.coord, groups=groups)
        if spec.nspec > 1:
            group_utils.show_group_meta(meta, idkey=sdb.idkey, show_all_keys=False)
            raise ValueError("Retreived more than 1 spectrum.  Choose your GROUP with --group=")
        spec_file = pargs.spec_file+'_{:s}'.format(meta['GROUP'][0])
    else:
        spec = pargs.spec_file
        spec_file = pargs.spec_file
    # Save spectrum filename to AbsSystem
    abs_sys.spec_file = spec_file

    # Run
    gui = XAbsSysGui(spec, abs_sys, norm=norm, llist=llist, outfil=pargs.outfile)
    gui.show()
    app.exec_()

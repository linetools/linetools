""" Module for XAbsSysGui
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import io
import json
import pdb

from PyQt4 import QtGui
from PyQt4 import QtCore

import warnings

from astropy import units as u
from linetools.isgm import utils as ltiu
from linetools.guis import line_widgets as ltgl

'''
=======
Analyzing system for future kinematic analysis

Here is now my preferred approach to perform the
analysis:

1.  Inspect the velocity plots.
2.  Identify the best low-ion transition for the analysis.
  a. High S/N
  b. Strong, but not saturated (or just barely)
  c. Preferably SiII, SII, ZnII, MgII (i.e. highly not refractory)
3.  Hit "^" on the line for a low-ion kinematic tracer
  a.  Adjust velocity limits if need be (1, 2)
4.  Hit "&" on the line for a high-ion kinematic tracer

See also VelPlotWidget doc

=======
Analyzing system for future abundance analysis

See VelPlotWidget doc
'''

class XAbsSysGui(QtGui.QDialog):
    """ GUI to replace XIDL x_velplot (and more)
    """
    def __init__(self, ispec, abs_sys, parent=None, llist=None, norm=True,
                 vmnx=[-300., 300.]*u.km/u.s, outfil=None):
        """
        spec : Filename or Spectrum1D
        abs_sys : AbsSystem
          Absorption system class
        Norm : bool, optional
          Normalized spectrum?
        """
        from linetools.guis import spec_widgets as ltgs
        super(XAbsSysGui, self).__init__(parent)

        # Initialize
        self.abs_sys = abs_sys
        self.z = self.abs_sys.zabs
        abs_lines = abs_sys.list_of_abslines()
        self.vmnx = vmnx
        if outfil is None:
            self.outfil = 'tmp_abskin.json'
            warnings.warn("Outfil not specified.  Using {:s} as the default".format(self.outfil))
        else:
            self.outfil = outfil
        self.norm = norm

        # Grab the pieces and tie together
        newfont = QtGui.QFont("Times", 10, QtGui.QFont.Bold)
        sys_label = QtGui.QLabel('Name: \n {:s}'.format(abs_sys.name))
        sys_label.setFont(newfont)
        self.vplt_widg = ltgs.VelPlotWidget(ispec, self.z, abs_lines=abs_lines, llist=llist,
                                            vmnx=self.vmnx, norm=self.norm)
        self.pltline_widg = ltgl.PlotLinesWidget(init_llist=self.vplt_widg.llist,
                                                 init_z=self.z, edit_z=False)
        #self.pltline_widg.spec_widg = self.vplt_widg

        self.slines = ltgl.SelectedLinesWidget(self.vplt_widg.llist[self.vplt_widg.llist['List']],
                                               init_select=self.vplt_widg.llist['show_line'],
                                               plot_widget=self.vplt_widg)

        # Connections
        self.pltline_widg.llist_widget.currentItemChanged.connect(self.on_llist_change)
        #self.connect(self.pltline_widg.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)
        self.vplt_widg.canvas.mpl_connect('key_press_event', self.on_key)

        # Outfil
        wbtn = QtGui.QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.write_out)
        self.out_box = QtGui.QLineEdit()
        self.out_box.setText(self.outfil)
        self.connect(self.out_box, QtCore.SIGNAL('editingFinished ()'), self.set_outfil)

        #QtCore.pyqtRemoveInputHook()
        #pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Quit
        buttons = QtGui.QWidget()
        wqbtn = QtGui.QPushButton('Write+Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.write_quit)
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.quit)

        # Sizes
        lines_widg = QtGui.QWidget()
        lines_widg.setMaximumWidth(300)
        lines_widg.setMinimumWidth(200)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(sys_label)
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(self.slines)
        vbox.addWidget(wbtn)
        vbox.addWidget(self.out_box)
        # Write/Quit buttons
        hbox1 = QtGui.QHBoxLayout()
        hbox1.addWidget(wqbtn)
        hbox1.addWidget(qbtn)
        buttons.setLayout(hbox1)
        #
        vbox.addWidget(buttons)
        lines_widg.setLayout(vbox)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.vplt_widg)
        hbox.addWidget(lines_widg)

        self.setLayout(hbox)
        # Initial draw
        self.vplt_widg.on_draw()

    # Overload, as needed
    def on_key(self, event):
        pass

    # Change list of lines to choose from
    def on_llist_change(self):
        llist = self.pltline_widg.llist
        all_lines = list( llist[llist['List']]._data['wrest'] )
        # Set selected
        wrest = [line.wrest for line in self.vplt_widg.abs_lines]
        select = []
        for iwrest in wrest:
            try:
                select.append(all_lines.index(iwrest))
            except ValueError:
                pass
        select.sort()
        # GUIs
        self.vplt_widg.llist['List'] = llist['List']
        self.vplt_widg.llist['show_line'] = select
        self.vplt_widg.idx_line = 0
        self.slines.selected = select
        self.slines.on_list_change(llist[llist['List']])

    # Write
    def set_outfil(self):
        self.outfil = str(self.out_box.text())
        print('AbsKin: Will write to {:s}'.format(self.outfil))

    '''
    # Set z from pltline_widg
    def setz(self):
        self.vplt_widg.z = self.pltline_widg.llist['z']
        self.z = self.pltline_widg.llist['z']
        self.vplt_widg.on_draw()
    '''

    def set_new_comps(self):
        """ Generate new components and fill into abs_sys
        Ignores velocity limits when building
        """
        # Add spectrum filename, coord
        abs_lines = self.vplt_widg.abs_lines
        for line in abs_lines:
            line.analy['datafile'] = self.vplt_widg.spec_fil
            line.attrib['coord'] = self.abs_sys.coord
        # Components
        comps = ltiu.build_components_from_abslines(abs_lines, chk_vel=False)
        self.abs_sys._components = comps
        # Return
        return

    # Write
    def write_out(self):
        # Update components and spectrum filename
        self.set_new_comps()
        # Dict
        adict = self.abs_sys.to_dict()

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        print("Wrote abs_sys to {:s}".format(self.outfil))
        with io.open(self.outfil, 'w', encoding='utf-8') as f:
            f.write(json.dumps(adict, sort_keys=True, indent=4,
                                       separators=(',', ': ')))

    # Write + Quit
    def write_quit(self):
        self.write_out()
        self.flg_quit = 1
        self.done(1)

    # Write + Quit
    def quit(self):
        self.flg_quit = 0
        self.done(1)


def main(ispec, abs_sys, **kwargs):
    """ Simple GUI call
    ParametersS
    ----------
    args
    kwargs

    Returns
    -------

    """
    import sys
    # Run
    app = QtGui.QApplication(sys.argv)
    gui = XAbsSysGui(ispec, abs_sys, **kwargs)
    gui.show()
    app.exec_()
    return

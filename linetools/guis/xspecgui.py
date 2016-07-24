""" Module for XSpecGui
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import sys
from PyQt4 import QtGui

from matplotlib import rcParams

from astropy.units import Quantity

from linetools.guis import utils as ltgu
from linetools.guis import line_widgets as ltgl
from linetools.guis import spec_widgets as ltgsp

class XSpecGui(QtGui.QMainWindow):
    """ GUI to replace XIDL x_specplot (which simulated a GUI by T. Barlow)
    """
    def __init__(self, ispec, parent=None, zsys=None, norm=None, exten=None,
                 rsp_kwargs={}, unit_test=False, **kwargs):
        QtGui.QMainWindow.__init__(self, parent)
        """
        ispec = str, XSpectrum1D or tuple of arrays
          Input spectrum or spectrum filename.  If tuple then (wave,
          fx), (wave, fx, sig) or (wave, fx, sig, co)
        parent : Widget parent, optional
        zsys : float, optional
          intial redshift
        exten : int, optional
          extension for the spectrum in multi-extension FITS file
        norm : bool, optional
          True if the spectrum is normalized
        """
        #reload(ltgl)
        #reload(ltgsp)
        # INIT
        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Needed to avoid crash in large spectral files
        rcParams['agg.path.chunksize'] = 20000
        rcParams['axes.formatter.useoffset'] = False  # avoid scientific notation in axes tick labels

        # Build a widget combining several others
        self.main_widget = QtGui.QWidget()

        # Status bar
        self.create_status_bar()

        # Grab the pieces and tie together
        self.pltline_widg = ltgl.PlotLinesWidget(status=self.statusBar, init_z=zsys)
        self.pltline_widg.setMaximumWidth(300)

        # Hook the spec widget to Plot Line
        self.spec_widg = ltgsp.ExamineSpecWidget(ispec,status=self.statusBar,
                                                 parent=self,
                                                llist=self.pltline_widg.llist,
                                                zsys=zsys, norm=norm, exten=exten,
                                                 rsp_kwargs=rsp_kwargs, **kwargs)
        self.pltline_widg.spec_widg = self.spec_widg

        self.spec_widg.canvas.mpl_connect('button_press_event', self.on_click)

        extras = QtGui.QWidget()
        extras.setMaximumWidth(130)
        vbox = QtGui.QVBoxLayout()
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.clicked.connect(self.quit)
        vbox.addWidget(self.pltline_widg)
        vbox.addWidget(qbtn)
        extras.setLayout(vbox)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(self.spec_widg)
        hbox.addWidget(extras)

        self.main_widget.setLayout(hbox)

        # Point MainWindow
        self.setCentralWidget(self.main_widget)
        if unit_test:
            self.quit()

    def create_status_bar(self):
        """ Status bar for the GUI
        """
        self.status_text = QtGui.QLabel("XSpec")
        self.statusBar().addWidget(self.status_text, 1)

    def on_click(self, event):
        """ Over-loads click events
        """
        if event.button == 3: # Set redshift
            if self.pltline_widg.llist['List'] is None:
                return
            self.select_line_widg = ltgl.SelectLineWidget(
                self.pltline_widg.llist[self.pltline_widg.llist['List']]._data)
            self.select_line_widg.exec_()
            line = self.select_line_widg.line
            if line.strip() == 'None':
                return
            #
            quant = line.split('::')[1].lstrip()
            spltw = quant.split(' ')
            wrest = Quantity(float(spltw[0]), unit=spltw[1])
            z = event.xdata/wrest.value - 1.
            self.pltline_widg.llist['z'] = z
            self.statusBar().showMessage('z = {:f}'.format(z))

            self.pltline_widg.zbox.setText('{:.5f}'.format(self.pltline_widg.llist['z']))

            # Draw
            self.spec_widg.on_draw()
    # Quit
    def quit(self):
        self.close()


def main(args, **kwargs):
    from linetools.spectra.xspectrum1d import XSpectrum1D

    if not isinstance(args,(XSpectrum1D,tuple,basestring)):
        raise IOError("Bad input")
    # Run
    app = QtGui.QApplication(sys.argv)
    gui = XSpecGui(args, **kwargs)
    gui.show()
    app.exec_()
    return

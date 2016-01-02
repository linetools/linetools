""" Module for line widets
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt4 import QtGui
from PyQt4 import QtCore

from astropy.table import Table

from linetools.guis import utils as ltgu

# #####
class PlotLinesWidget(QtGui.QWidget):
    """ Widget to set up spectral lines for plotting
    """
    def __init__(self, parent=None, status=None, init_llist=None, init_z=None):
        """
        Parameters
        ----------
        parent
        status
        init_llist
        init_z

        Returns
        -------

        """
        super(PlotLinesWidget, self).__init__(parent)

        # Initialize
        if not status is None:
            self.statusBar = status
        if init_z is None:
            init_z = 0.

        # Create a dialog window for redshift
        z_label = QtGui.QLabel('z=')
        self.zbox = QtGui.QLineEdit()
        self.zbox.z_frmt = '{:.7f}'
        self.zbox.setText(self.zbox.z_frmt.format(init_z))
        self.zbox.setMinimumWidth(50)
        self.connect(self.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)

        # Create the line list
        self.lists = ['None', 'ISM', 'Strong', 'Galaxy', 'H2', 'EUV', 'OVI']
        #'grb.lst', 'dla.lst', 'lls.lst', 'subLLS.lst',
#                      'lyman.lst', 'Dlyman.lst', 'gal_vac.lst', 'ne8.lst',
#                      'lowz_ovi.lst', 'casbah.lst', 'H2.lst']
        list_label = QtGui.QLabel('Line Lists:')
        self.llist_widget = QtGui.QListWidget(self)
        for ilist in self.lists:
            self.llist_widget.addItem(ilist)
        self.llist_widget.setCurrentRow(0)
        self.llist_widget.currentItemChanged.connect(self.on_list_change)
        self.llist_widget.setMaximumHeight(100)

        # Input line list?
        if init_llist is None:
            self.llist = {} # Dict for the line lists
            self.llist['Plot'] = False
            self.llist['z'] = 0.
            self.llist['List'] = 'None'
            self.llist['Lists'] = []
        else: # Fill it all up and select
            self.llist = init_llist
            if not init_llist['List'] in self.lists:
                self.lists.append(init_llist['List'])
                self.llist_widget.addItem(init_llist['List'])
                self.llist_widget.setCurrentRow(len(self.lists)-1)
            else:
                idx = self.lists.index(init_llist['List'])
                self.llist_widget.setCurrentRow(idx)
            try:
                self.zbox.setText(self.zbox.z_frmt.format(init_llist['z']))
            except KeyError:
                pass

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(self.zbox)
        vbox.addWidget(list_label)
        vbox.addWidget(self.llist_widget)

        self.setLayout(vbox)
        self.setMaximumHeight(200)

    def on_list_change(self,curr,prev):
        llist = str(curr.text())
        # Print
        try:
            self.statusBar().showMessage('You chose: {:s}'.format(llist))
        except AttributeError:
            print('You chose: {:s}'.format(curr.text()))

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()
        self.llist = ltgu.set_llist(llist,in_dict=self.llist)

        # Try to draw
        if self.llist['Plot'] is True:
            try:
                self.spec_widg.on_draw()
            except AttributeError:
                return

    def setz(self):
        sstr = unicode(self.zbox.text())
        try:
            self.llist['z'] = float(sstr)
        except ValueError:
            try:
                self.statusBar().showMessage('ERROR: z Input must be a float! Try again..')
            except AttributeError:
                print('ERROR: z Input must be a float! Try again..')
            self.zbox.setText(self.zbox.z_frmt.format(self.llist['z']))
            return

        # Report
        try:
            self.statusBar().showMessage('z = {:g}'.format(self.llist['z']))
        except AttributeError:
            print('z = {:g}'.format(self.llist['z']))

        # Try to draw
        try:
            self.spec_widg.on_draw()
        except AttributeError:
            return

class SelectLineWidget(QtGui.QDialog):
    """ Widget to select a spectral line
    inp: string or dict or Table
      Input line list

    15-Dec-2014 by JXP
    """
    def __init__(self, inp, parent=None):
        super(SelectLineWidget, self).__init__(parent)

        # Line list Table
        if isinstance(inp, Table):
            lines = inp
        else:
            raise ValueError('SelectLineWidget: Wrong type of input')

        self.resize(250, 800)

        # Create the line list
        line_label = QtGui.QLabel('Lines:')
        self.lines_widget = QtGui.QListWidget(self)
        self.lines_widget.addItem('None')
        self.lines_widget.setCurrentRow(0)

        #xdb.set_trace()
        # Loop on lines (could put a preferred list first)
        # Sort
        srt = np.argsort(lines['wrest'])
        for ii in srt:
            self.lines_widget.addItem('{:s} :: {:.3f} :: {}'.format(lines['name'][ii],
                                                         lines['wrest'][ii], lines['f'][ii]))
        self.lines_widget.currentItemChanged.connect(self.on_list_change)
        #self.scrollArea = QtGui.QScrollArea()

        # Quit
        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.clicked.connect(self.close)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(line_label)
        vbox.addWidget(self.lines_widget)
        vbox.addWidget(qbtn)

        self.setLayout(vbox)

    def on_list_change(self, curr, prev):
        self.line = str(curr.text())
        # Print
        print('You chose: {:s}'.format(curr.text()))
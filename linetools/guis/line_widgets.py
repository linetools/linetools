""" Module for line widets
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt4 import QtGui
from PyQt4 import QtCore

import numpy as np
import pdb

from astropy.table import Table

from linetools.guis import utils as ltgu
from linetools.lists.linelist import LineList


# #####
class PlotLinesWidget(QtGui.QWidget):
    """ Widget to set up spectral lines for plotting
    """
    def __init__(self, parent=None, status=None, init_llist=None, init_z=None,
                 edit_z=True):
        """
        Parameters
        ----------
        parent : Widget parent
        status : Point to status bar
        init_llist : input LineList dictionary (from another widget)
        init_z : float, optional
          Initial redshift
        edit_z : bool, optional
          Allow z to be editable

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
        if edit_z:
            z_label = QtGui.QLabel('z=')
            self.zbox = QtGui.QLineEdit()
            self.zbox.z_frmt = '{:.7f}'
            self.zbox.setText(self.zbox.z_frmt.format(init_z))
            self.zbox.setMinimumWidth(50)
            self.connect(self.zbox, QtCore.SIGNAL('editingFinished ()'), self.setz)
        else:
            z_label = QtGui.QLabel('z={:.7f}'.format(init_z))

        # Create the line list
        self.lists = ['None', 'ISM', 'Strong', 'HI', 'Galaxy', 'H2', 'EUV']
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
            self.llist['List'] = 'None'  # Name of the LineList being used
            self.llist['Lists'] = []     # Archived LineLists
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
            except (AttributeError, KeyError):
                pass

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        if edit_z:
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

        # Loop on lines (could put a preferred list first)
        # Sort
        srt = np.argsort(lines['wrest'])
        for ii in srt:
            try:
                s = '{:s} :: {:.2f} :: {:.3f}'.format(lines['name'][ii], lines['wrest'][ii],
                                                      lines['f'][ii])
            except ValueError:  # f-value masked (most likely)
                s = '{:s} :: {:.2f}'.format(lines['name'][ii], lines['wrest'][ii])
            #  is there a column called 'redshift'? (only used in igmguesses for now)
            try:
                s += ' :: z{:.3f}'.format(lines['redshift'][ii])
                self.resize(350, 800)
            except KeyError:
                pass
            self.lines_widget.addItem(s)

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


class SelectedLinesWidget(QtGui.QWidget):
    """ Widget to show and enable lines to be selected
    inp : LineList
      Input LineList
    init_select : str or list of indices
      str -- 'All'

    """
    def __init__(self, inp, parent=None, init_select=None, plot_widget=None):
        """
        """
        super(SelectedLinesWidget, self).__init__(parent)

        self.parent=parent

        # Line list Table
        if isinstance(inp, LineList):
            self.lines = inp._data
            self.llst = inp
        elif isinstance(inp,Table):
            raise ValueError('SelectedLineWidget: DEPRECATED')
        else:
            raise ValueError('SelectedLineWidget: Wrong type of input')

        self.plot_widget = plot_widget

        # Create the line list
        line_label = QtGui.QLabel('Lines:')
        self.lines_widget = QtGui.QListWidget(self)
        self.lines_widget.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)

        # Initialize list
        self.item_flg = 0
        self.init_list()

        # Initial selection
        if init_select is None:
            self.selected = [0]
        elif init_select == 'All':
            self.selected = []
            for ii in range(self.lines_widget.count()):
                self.lines_widget.item(ii).setSelected(True)
                self.selected.append(ii)
        else:
            self.selected = init_select
            if len(self.selected) == 0:
                self.selected = [0]

        for iselect in self.selected:
            self.lines_widget.item(iselect).setSelected(True)

        self.lines_widget.scrollToItem( self.lines_widget.item( self.selected[0] ) )

        # Events
        self.lines_widget.itemSelectionChanged.connect(self.on_item_change)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(line_label)
        vbox.addWidget(self.lines_widget)

        self.setLayout(vbox)

    def init_list(self):
        nlin = len(self.lines['wrest'])
        for ii in range(nlin):
            self.lines_widget.addItem('{:s} :: {:.3f} :: {}'.format(self.lines['name'][ii],
                                                         self.lines['wrest'][ii].value,
                                                         self.lines['f'][ii]))

    def on_item_change(self): #,item):
        # For big changes
        if self.item_flg == 1:
            return
        all_items = [self.lines_widget.item(ii) for ii in range(self.lines_widget.count())]
        sel_items = self.lines_widget.selectedItems()
        self.selected = [all_items.index(isel) for isel in sel_items]
        self.selected.sort()

        #QtCore.pyqtRemoveInputHook()
        #xdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # Update llist
        try:
            self.plot_widget.llist['show_line'] = self.selected
        except AttributeError:
            if self.parent is not None:
                self.parent.updated_slines(self.selected)
            return
        else:
            self.plot_widget.on_draw()

    def on_list_change(self, llist):
        # Clear
        if not isinstance(llist, LineList):
            raise ValueError('Expecting LineList!!')
        self.item_flg = 1
        self.lines = llist._data
        self.llst = llist
        self.lines_widget.clear()
        # Initialize
        self.init_list()
        # Set selected
        for iselect in self.selected:
            self.lines_widget.item(iselect).setSelected(True)
        self.lines_widget.scrollToItem(self.lines_widget.item(self.selected[0]))
        self.item_flg = 0

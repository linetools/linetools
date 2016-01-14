""" Module for simple widgets
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt4 import QtGui
from PyQt4 import QtCore


class AnsBox(QtGui.QDialog):
    """Solicit an input answer from the User
    """
    def __init__(self, lbl, format=str, parent=None):
        """
        Parameters
        ----------
        lbl : str
        format : str
          Format for value
        """
        super(AnsBox, self).__init__(parent)

        self.format=format
        #
        label = QtGui.QLabel(lbl)
        self.box = QtGui.QLineEdit()
        self.box.setMinimumWidth(90)
        # Connect
        self.connect(self.box,
            QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self):
        try:
            self.value = self.format(unicode(self.box.text()))
        except ValueError:
            print('Bad input value! Try again with right type')
        else:
            self.done(0)

class EditBox(QtGui.QWidget):
    """ Generate a simple box for editing

    """
    def __init__(self, initv, lbl, format, parent=None):
        """
        Parameters
        ----------
        initv : str
          Initial value
        lbl : str
        format : str
          Format for text
        """
        super(EditBox, self).__init__(parent)

        self.value = initv
        #
        label = QtGui.QLabel(lbl)
        self.box = QtGui.QLineEdit()
        # Format
        self.box.frmt = format
        self.box.setText(self.box.frmt.format(self.value))
        self.box.setMinimumWidth(90)
        # Connect
        self.connect(self.box,
            QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self):
        self.value = unicode(self.box.text())

    def set_text(self,value):
        self.value = value
        self.box.setText(self.box.frmt.format(self.value))


# ##################################
class WarningWidg(QtGui.QDialog):
    """ GUI to warn user about coming action and solicit response
    """
    def __init__(self, message, parent=None):
        """
        Parameters
        ----------
        message : str
          Message to display
        """
        super(WarningWidg, self).__init__(parent)

        # Initialize

        # Grab the pieces and tie together
        z_label = QtGui.QLabel('Warning: {:s}'.format(message))

        # Quit
        nbtn = QtGui.QPushButton('No', self)
        nbtn.clicked.connect(self.touch_no)
        ybtn = QtGui.QPushButton('Yes', self)
        ybtn.clicked.connect(self.touch_yes)

        # Layout
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(z_label)
        vbox.addWidget(nbtn)
        vbox.addWidget(ybtn)
        self.setLayout(vbox)

    def touch_yes(self):
        self.ans = True
        self.done(0)

    def touch_no(self):
        self.ans = False
        self.done(0)


class WriteQuitWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        """
        """
        super(WriteQuitWidget, self).__init__(parent)
        self.parent = parent

        # Generate Buttons
        wbtn = QtGui.QPushButton('Write', self)
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.parent.write_out)

        wqbtn = QtGui.QPushButton('Write\n Quit', self)
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.parent.write_quit)

        qbtn = QtGui.QPushButton('Quit', self)
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.parent.quit)

        # Layout
        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(wbtn)
        hbox.addWidget(wqbtn)
        hbox.addWidget(qbtn)
        self.setLayout(hbox)


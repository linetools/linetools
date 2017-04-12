""" Module for simple widgets
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtWidgets import QWidget, QDialog, QLabel
from PyQt5.QtWidgets import QPushButton, QLineEdit
from PyQt5.QtWidgets import QVBoxLayout, QHBoxLayout

try:
    ustr = unicode
except NameError:  # For Python 3
    ustr = str

class AnsBox(QDialog):
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
        label = QLabel(lbl)
        self.box = QLineEdit()
        self.box.setMinimumWidth(90)
        # Connect
        self.box.textChanged[str].connect(self.setv)
        self.box.editingFinished.connect(self.finish)
        # Layout
        vbox = QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self, text):
        self.box.setText(text)
        self.box.adjustSize()

    def finish(self):
        try:
            self.value = self.format(ustr(self.box.text()))
        except ValueError:
            print('Bad input value! Try again with right type')
        else:
            self.done(0)

class EditBox(QWidget):
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
        label = QLabel(lbl)
        self.box = QLineEdit()
        # Format
        self.box.frmt = format
        self.box.setText(self.box.frmt.format(self.value))
        self.box.setMinimumWidth(90)
        # Connect
        self.box.textChanged[str].connect(self.setv)
        #self.connect(self.box,
        #    QtCore.SIGNAL('editingFinished ()'), self.setv)
        # Layout
        vbox = QVBoxLayout()
        vbox.addWidget(label)
        vbox.addWidget(self.box)
        self.setLayout(vbox)

    def setv(self, text):
        self.box.setText(text)
        self.box.adjustSize()
        self.value = str(self.box.text())

    def set_text(self,value):
        self.value = value
        self.box.setText(self.box.frmt.format(self.value))


# ##################################
class WarningWidg(QDialog):
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
        z_label = QLabel('Warning: {:s}'.format(message))

        # Quit
        nbtn = QPushButton(self)
        nbtn.setText('No')
        nbtn.clicked.connect(self.touch_no)
        ybtn = QPushButton(self)
        ybtn.setText('Yes')
        ybtn.clicked.connect(self.touch_yes)

        # Layout
        vbox = QVBoxLayout()
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


class WriteQuitWidget(QWidget):
    def __init__(self, parent=None):
        """
        """
        super(WriteQuitWidget, self).__init__(parent)
        self.parent = parent

        # Generate Buttons
        wbtn = QPushButton(self)
        wbtn.setText('Write')
        wbtn.setAutoDefault(False)
        wbtn.clicked.connect(self.parent.write_out)

        wqbtn = QPushButton(self)
        wqbtn.setText('Write\n Quit')
        wqbtn.setAutoDefault(False)
        wqbtn.clicked.connect(self.parent.write_quit)

        qbtn = QPushButton(self)
        qbtn.setText('Quit')
        qbtn.setAutoDefault(False)
        qbtn.clicked.connect(self.parent.quit)

        # Layout
        hbox = QHBoxLayout()
        hbox.addWidget(wbtn)
        hbox.addWidget(wqbtn)
        hbox.addWidget(qbtn)
        self.setLayout(hbox)


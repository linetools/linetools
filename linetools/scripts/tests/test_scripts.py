# Module to test scripts from linetools
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS
import imp

import pytest
from linetools.scripts.lt_absline import plot_absline
from linetools.scripts import lt_line
from linetools.scripts import lt_plot

def test_lt_absline():
    plot_absline(1550, 15, 5, show=False)
    plot_absline('CIV 1548', 15, 5, show=False)


def test_lt_line():
    # Run through the motions
    lt_line.main(['HI'])
    lt_line.main(['HI1215'])
    lt_line.main(['1215'])
    lt_line.main(['--all'])
    #lt_line.main()
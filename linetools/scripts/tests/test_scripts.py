# Module to test scripts from linetools
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pdb
from linetools.scripts.lt_absline import plot_absline
from linetools.scripts import lt_line


def test_lt_absline():
    plot_absline(1550, 15, 5, show=False)

def test_lt_line():
    # Run through the motions
    lt_line.main(['HI'])
    lt_line.main(['HI1215'])
    lt_line.main(['1215'])
    lt_line.main(['--all'])

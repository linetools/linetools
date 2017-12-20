# Module to test scripts from linetools
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
from linetools.scripts.lt_absline import plot_absline
from linetools.scripts import lt_line
from linetools.scripts import lt_radec
from linetools.scripts import lt_solabnd
from linetools.scripts import lt_get_COS_LP


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


def test_lt_radec():
    lt_radec.main(['152.25900,7.22885'])
    lt_radec.main(['J100902.16+071343.8'])
    lt_radec.main(['10:09:02.16,+07:13:43.8'])

def test_lt_solabnd():
    lt_solabnd.main(['Fe'])
    lt_solabnd.main(['-a'])
    lt_solabnd.main(['-a', '--sortZ'])

def test_lt_get_COS_LP():
    lt_get_COS_LP.main(["2017-10-01"])

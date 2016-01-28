from __future__ import print_function, absolute_import, \
     division, unicode_literals

from linetools.scripts.lt_absline import plot_absline


def test_lt_absline():
    plot_absline(1550, 15, 5, show=False)


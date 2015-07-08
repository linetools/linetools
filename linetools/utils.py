""" Module for general utilities which don't belong an a sub-package.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

def between(a, vmin, vmax):
    """ Return a boolean array True where vmin <= a < vmax.

    Notes
    -----
    Careful of floating point issues when dealing with equalities.
    """
    a = np.asarray(a)
    c = a < vmax
    c &= a >= vmin
    return c

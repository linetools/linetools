""" General utilities.
"""

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

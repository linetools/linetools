""" Module for general utilities which don't belong in a sub-package.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

def between(a, vmin, vmax):
    """ Return a boolean array True where vmin <= a < vmax.

    Parameters
    ----------
    a : array of shape N
      Input array.
    vmin, vmix: float
      Minimum and maximum values to test between.

    Returns
    -------
    c : array of shape N
      Boolean array true where vmin < a < vmax.

    Notes
    -----
    This is a convenience function equivalent to (vmin <= a) & (a < vmax).
    Be careful of floating point issues when dealing with equalities.
    """
    a = np.asarray(a)
    c = a < vmax
    c &= a >= vmin
    return c

#
def scipy_rebin(a, *args):
    ''' Simple script to rebin an input array to a new shape.  Akin to IDL's routine
    Taken from scipy documentation: http://wiki.scipy.org/Cookbook/Rebinning
    As in IDL, the new shape must be a factor of the old one.
    The ugly 'evList trick' builds and executes a python command 
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))


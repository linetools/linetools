""" Plotting utilities."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

def get_flux_plotrange(fl):
    """ Find y limits for a spectrum. 
    """
    ymax = abs(np.percentile(fl[~np.isnan(fl)], 95)) * 1.5
    return -0.1 * ymax, ymax

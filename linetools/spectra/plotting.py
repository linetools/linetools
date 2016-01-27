""" Plotting utilities."""
from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np
from astropy.units.quantity import Quantity

def get_flux_plotrange(fl, perc=95, mult_pos=1.5, mult_neg=-0.1):
    """ Find y limits for a spectrum. 
    """
    if isinstance(fl, Quantity):
        fl = fl.value
    else:
        fl = np.array(fl, copy=True)

    try: 
        ymax = abs(np.percentile(fl[~np.isnan(fl)], perc)) * mult_pos
    except:
        import pdb; pdb.set_trace()
    return mult_neg * ymax, ymax

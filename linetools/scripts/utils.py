""" Utils for scripts
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb

def coord_arg_to_coord(carg):
    """
    Parameters
    ----------
    carg : str
      Argument from parser for coordinates
      Eligible formats are like:
         J081240.7+320809
         122.223,-23.2322
         07:45:00.47,34:17:31.1

    Returns
    -------
    icoord : str or tuple
      Allowed format for coord input to linetools.utils.radec_to_coord

    """
    if ',' in carg:
        radec = carg.split(',')
        if ':' in radec[0]:   # 07:45:00.47,34:17:31.1
            icoord = (radec[0], radec[1])
        else:  # 122.223,-23.2322
            icoord = (float(radec[0]), float(radec[1]))
    else:  # J081240.7+320809
        icoord = carg
    # Return
    return icoord


def get_today_str():
    """Returns a string representation of current day
    format YYYY-MM-DD"""
    import time
    t = time.gmtime()
    year = str(t.tm_year)
    month = str(t.tm_mon)
    day = str(t.tm_mday)
    if len(month) == 1:
        month = '0' + month
    if len(day) == 1:
        day = '0' + day
    s = '{}-{}-{}'.format(year, month, day)
    return s
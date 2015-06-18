"""
Module for LineList Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy import constants as const
from astropy.io import fits

#from xastropy.xutils import xdebug as xdb

#
class LineList(Spectrum1D):
    '''Class to over-load Spectrum1D for new functionality not yet in specutils
    '''


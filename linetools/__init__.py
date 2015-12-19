# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, absolute_import, division, unicode_literals


"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    pass
    #from example_mod import *
    #from . import spectra
    #from . import lists

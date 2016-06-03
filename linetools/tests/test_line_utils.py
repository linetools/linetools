# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, pdb
from astropy import units as u

from linetools.spectralline import AbsLine
from linetools import line_utils as ltlu


def test_parse_abslines():
    # Init AbsLines
    abslines = [AbsLine(1215.6700*u.AA), AbsLine('CII 1334')]
    # wrest
    wrest_values = ltlu.parse_abslines(abslines, 'wrest')
    np.testing.assert_allclose(wrest_values[1], 1334.5323*u.AA)
    # EW
    EW_values = ltlu.parse_abslines(abslines, 'EW')
    np.testing.assert_allclose(EW_values[1].value, 0.)
    # data
    A_values = ltlu.parse_abslines(abslines, 'A')
    np.testing.assert_allclose(A_values[0].value, 626500000.0)

def test_transtabl():
    # Init AbsLines
    abslines = [AbsLine(1215.6700*u.AA), AbsLine('CII 1334')]
    #
    tbl = ltlu.transtable_from_abslines(abslines)
    assert len(tbl) == 2

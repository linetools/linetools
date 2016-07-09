# Module to test scripts from linetools
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest
from linetools.scripts.lt_absline import plot_absline
from linetools.scripts import lt_line, lt_viewhdf
from linetools.spectra.xspectrum1d import XSpectrum1D


def test_lt_absline():
    plot_absline(1550, 15, 5, show=False)
    plot_absline('CIV 1548', 15, 5, show=False)


def test_lt_line():
    # Run through the motions
    lt_line.main(['HI'])
    lt_line.main(['HI1215'])
    lt_line.main(['1215'])
    lt_line.main(['--all'])


def test_viewhdf():
    # Generate an HDF5
    wv = np.arange(1000)
    fx = np.ones_like(wv)
    sig = np.ones_like(wv)
    spec = XSpectrum1D.from_tuple((wv,fx,sig))
    spec.write_to_hdf5('tmp.hdf5')
    # Tickle
    lt_viewhdf.main(['tmp.hdf5','-f'])
    # Cleanup
    os.remove('tmp.hdf5')

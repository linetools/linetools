#!/usr/bin/env python

"""
Plot an absorption line with given parameters.
Also print the line data
  Examples:
  linet_absline.py 1334.3    13.5 7
  linet_absline.py 1215.6701 14.0 30
"""

import sys
import os
import numpy as np
import pdb
import argparse

from astropy import units as u

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.lists.linelist import LineList
from linetools.spectralline import AbsLine
from linetools.analysis import voigt as lav

# Parse
parser = argparse.ArgumentParser(description='Parser for linet_absline')
parser.add_argument("wrest", type=float, help="Rest wavelength in Angstroms")
parser.add_argument("logN", type=float, help="log10 column density")
parser.add_argument("b", type=float, help="b-value in km/s")
 
pargs = parser.parse_args()

# Search for the closest absline
aline = AbsLine(pargs.wrest*u.AA, closest=True)

# Generate a fake wavelength array near the line
wvoff = 50. # Ang
dwv = pargs.wrest/100000. # Ang (echelle)
wave = np.arange(pargs.wrest-wvoff, pargs.wrest+wvoff, dwv)

# Generate spectrum with voigt
aline.attrib['N'] = pargs.logN
aline.attrib['b'] = pargs.b * u.km/u.s
xspec = aline.generate_voigt(wave=wave*u.AA)
xspec.dummy_sig(0.1) # S/N = 10 per pix


# Calculate EW
aline.analy['spec'] = xspec
aline.analy['wvlim'] = np.array([pargs.wrest-15., pargs.wrest+15])*u.AA
aline.measure_ew()
print(aline)
print('EW = {:g}'.format(aline.attrib['EW']))

# Plot
xspec.plot()